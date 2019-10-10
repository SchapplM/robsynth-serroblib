% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPPRR1
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:25
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPPRR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:52
	% EndTime: 2019-10-09 23:25:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:52
	% EndTime: 2019-10-09 23:25:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:52
	% EndTime: 2019-10-09 23:25:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t32 = qJ(1) + pkin(9);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [(-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(1), 0, 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:52
	% EndTime: 2019-10-09 23:25:52
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (22->9), mult. (24->12), div. (0->0), fcn. (14->4), ass. (0->6)
	t13 = -pkin(2) + r_i_i_C(2);
	t12 = r_i_i_C(3) + qJ(3);
	t11 = qJ(1) + pkin(9);
	t10 = cos(t11);
	t9 = sin(t11);
	t1 = [t10 * qJD(3) + (-cos(qJ(1)) * pkin(1) - t12 * t9 + t13 * t10) * qJD(1), 0, qJD(1) * t10, 0, 0, 0; t9 * qJD(3) + (-sin(qJ(1)) * pkin(1) + t13 * t9 + t12 * t10) * qJD(1), 0, qJD(1) * t9, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:52
	% EndTime: 2019-10-09 23:25:52
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (33->12), mult. (32->14), div. (0->0), fcn. (20->4), ass. (0->8)
	t15 = r_i_i_C(2) + qJ(3);
	t12 = qJ(1) + pkin(9);
	t10 = sin(t12);
	t14 = qJD(1) * t10;
	t13 = -pkin(2) - r_i_i_C(3) - qJ(4);
	t11 = cos(t12);
	t9 = qJD(1) * t11;
	t1 = [t11 * qJD(3) - t10 * qJD(4) + (-cos(qJ(1)) * pkin(1) - t15 * t10 + t13 * t11) * qJD(1), 0, t9, -t14, 0, 0; t10 * qJD(3) + t11 * qJD(4) + (-sin(qJ(1)) * pkin(1) + t15 * t11 + t13 * t10) * qJD(1), 0, t14, t9, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:52
	% EndTime: 2019-10-09 23:25:53
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (64->23), mult. (84->35), div. (0->0), fcn. (54->6), ass. (0->16)
	t25 = sin(qJ(5));
	t26 = cos(qJ(5));
	t28 = (r_i_i_C(1) * t26 - r_i_i_C(2) * t25) * qJD(5);
	t36 = qJD(4) + t28;
	t24 = qJ(1) + pkin(9);
	t22 = sin(t24);
	t35 = qJD(1) * t22;
	t34 = qJD(1) * t25;
	t33 = qJD(1) * t26;
	t32 = qJD(5) * t25;
	t31 = qJD(5) * t26;
	t30 = pkin(7) + r_i_i_C(3) - qJ(3);
	t27 = -r_i_i_C(1) * t25 - r_i_i_C(2) * t26 - pkin(2) - qJ(4);
	t23 = cos(t24);
	t21 = qJD(1) * t23;
	t1 = [t23 * qJD(3) - t36 * t22 + (-cos(qJ(1)) * pkin(1) + t30 * t22 + t27 * t23) * qJD(1), 0, t21, -t35, (t22 * t34 - t23 * t31) * r_i_i_C(2) + (-t22 * t33 - t23 * t32) * r_i_i_C(1), 0; t22 * qJD(3) + t36 * t23 + (-sin(qJ(1)) * pkin(1) - t30 * t23 + t27 * t22) * qJD(1), 0, t35, t21, (-t22 * t31 - t23 * t34) * r_i_i_C(2) + (-t22 * t32 + t23 * t33) * r_i_i_C(1), 0; 0, 0, 0, 0, -t28, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:54
	% EndTime: 2019-10-09 23:25:54
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (186->43), mult. (290->66), div. (0->0), fcn. (225->8), ass. (0->34)
	t206 = cos(qJ(5));
	t204 = sin(qJ(5));
	t227 = pkin(8) + r_i_i_C(3);
	t229 = t227 * t204;
	t234 = (pkin(5) * t206 + t229) * qJD(5) + qJD(4);
	t203 = sin(qJ(6));
	t205 = cos(qJ(6));
	t212 = r_i_i_C(1) * t205 - r_i_i_C(2) * t203 + pkin(5);
	t232 = t212 * t206 + t229;
	t214 = qJD(1) * t204 + qJD(6);
	t231 = t205 * t214;
	t221 = qJD(6) * t204;
	t215 = qJD(1) + t221;
	t222 = qJD(5) * t206;
	t228 = t203 * t222 + t215 * t205;
	t225 = pkin(7) - qJ(3);
	t202 = qJ(1) + pkin(9);
	t200 = sin(t202);
	t224 = qJD(1) * t200;
	t201 = cos(t202);
	t199 = qJD(1) * t201;
	t223 = qJD(5) * t204;
	t220 = qJD(6) * t206;
	t217 = t227 * t206;
	t213 = r_i_i_C(1) * t203 + r_i_i_C(2) * t205;
	t211 = t214 * t203;
	t210 = -pkin(5) * t204 - pkin(2) - qJ(4) + t217;
	t209 = t215 * t203 - t205 * t222;
	t207 = -t213 * t220 + (-t212 * t204 + t217) * qJD(5);
	t198 = t209 * t200 - t201 * t231;
	t197 = t228 * t200 + t201 * t211;
	t196 = t200 * t231 + t209 * t201;
	t195 = t200 * t211 - t228 * t201;
	t1 = [t198 * r_i_i_C(1) + t197 * r_i_i_C(2) + t201 * qJD(3) - t234 * t200 + (-cos(qJ(1)) * pkin(1) + t225 * t200 + t210 * t201) * qJD(1), 0, t199, -t224, t207 * t201 - t224 * t232, t195 * r_i_i_C(1) + t196 * r_i_i_C(2); -t196 * r_i_i_C(1) + t195 * r_i_i_C(2) + t200 * qJD(3) + t234 * t201 + (-sin(qJ(1)) * pkin(1) - t225 * t201 + t210 * t200) * qJD(1), 0, t224, t199, t232 * t199 + t207 * t200, -t197 * r_i_i_C(1) + t198 * r_i_i_C(2); 0, 0, 0, 0, -qJD(5) * t232 + t213 * t221, (t203 * t220 + t205 * t223) * r_i_i_C(2) + (t203 * t223 - t205 * t220) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end