% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPRR12
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPPRR12_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPRR12_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR12_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:07:49
	% EndTime: 2019-12-31 18:07:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:07:49
	% EndTime: 2019-12-31 18:07:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:07:49
	% EndTime: 2019-12-31 18:07:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->6), mult. (20->10), div. (0->0), fcn. (12->2), ass. (0->5)
	t10 = -pkin(1) + r_i_i_C(2);
	t9 = r_i_i_C(3) + qJ(2);
	t8 = cos(qJ(1));
	t7 = sin(qJ(1));
	t1 = [t8 * qJD(2) + (t10 * t8 - t9 * t7) * qJD(1), qJD(1) * t8, 0, 0, 0; t7 * qJD(2) + (t10 * t7 + t9 * t8) * qJD(1), qJD(1) * t7, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:07:49
	% EndTime: 2019-12-31 18:07:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (15->10), mult. (36->14), div. (0->0), fcn. (24->4), ass. (0->7)
	t14 = sin(qJ(1));
	t18 = qJD(1) * t14;
	t17 = -pkin(1) - r_i_i_C(3) - qJ(3);
	t16 = r_i_i_C(1) * sin(pkin(8)) + r_i_i_C(2) * cos(pkin(8)) + qJ(2);
	t15 = cos(qJ(1));
	t11 = qJD(1) * t15;
	t1 = [t15 * qJD(2) - t14 * qJD(3) + (-t16 * t14 + t17 * t15) * qJD(1), t11, -t18, 0, 0; t14 * qJD(2) + t15 * qJD(3) + (t17 * t14 + t16 * t15) * qJD(1), t18, t11, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:07:49
	% EndTime: 2019-12-31 18:07:49
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (48->21), mult. (82->32), div. (0->0), fcn. (54->5), ass. (0->14)
	t24 = pkin(8) + qJ(4);
	t21 = sin(t24);
	t22 = cos(t24);
	t36 = (r_i_i_C(1) * t22 - r_i_i_C(2) * t21) * qJD(4);
	t27 = sin(qJ(1));
	t35 = qJD(1) * t27;
	t28 = cos(qJ(1));
	t23 = qJD(1) * t28;
	t34 = qJD(4) * t27;
	t33 = qJD(4) * t28;
	t32 = -pkin(1) - r_i_i_C(3) - pkin(6) - qJ(3);
	t30 = pkin(3) * sin(pkin(8)) + r_i_i_C(1) * t21 + r_i_i_C(2) * t22 + qJ(2);
	t29 = qJD(2) + t36;
	t1 = [-t27 * qJD(3) + t29 * t28 + (-t27 * t30 + t28 * t32) * qJD(1), t23, -t35, (-t21 * t23 - t22 * t34) * r_i_i_C(2) + (-t21 * t34 + t22 * t23) * r_i_i_C(1), 0; t28 * qJD(3) + t29 * t27 + (t27 * t32 + t28 * t30) * qJD(1), t35, t23, (-t21 * t35 + t22 * t33) * r_i_i_C(2) + (t21 * t33 + t22 * t35) * r_i_i_C(1), 0; 0, 0, 0, -t36, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:07:50
	% EndTime: 2019-12-31 18:07:50
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (165->45), mult. (288->78), div. (0->0), fcn. (225->7), ass. (0->33)
	t211 = sin(qJ(5));
	t213 = cos(qJ(5));
	t217 = (r_i_i_C(1) * t211 + r_i_i_C(2) * t213) * qJD(5);
	t232 = pkin(7) + r_i_i_C(3);
	t238 = t232 * qJD(4) - t217;
	t214 = cos(qJ(1));
	t208 = pkin(8) + qJ(4);
	t205 = sin(t208);
	t222 = qJD(1) * t205 + qJD(5);
	t237 = t222 * t214;
	t236 = t232 * t205;
	t206 = cos(t208);
	t234 = t232 * t206 - pkin(3) * sin(pkin(8)) - pkin(4) * t205 - qJ(2);
	t212 = sin(qJ(1));
	t227 = qJD(4) * t214;
	t233 = -t206 * t227 + t222 * t212;
	t231 = -pkin(1) - pkin(6) - qJ(3);
	t230 = qJD(1) * t212;
	t229 = qJD(4) * t212;
	t228 = qJD(4) * t213;
	t226 = qJD(5) * t206;
	t225 = qJD(1) * t232;
	t223 = -qJD(5) * t205 - qJD(1);
	t220 = r_i_i_C(1) * t213 - r_i_i_C(2) * t211 + pkin(4);
	t219 = t223 * t214;
	t216 = qJD(1) * t220;
	t215 = qJD(2) + (pkin(4) * t206 + t236) * qJD(4);
	t207 = qJD(1) * t214;
	t204 = t213 * t237 + (t206 * t228 + t223 * t211) * t212;
	t203 = t223 * t213 * t212 + (-t206 * t229 - t237) * t211;
	t202 = t211 * t219 - t233 * t213;
	t201 = t233 * t211 + t213 * t219;
	t1 = [t202 * r_i_i_C(1) + t201 * r_i_i_C(2) - t212 * qJD(3) + t215 * t214 + (t234 * t212 + t231 * t214) * qJD(1), t207, -t230, (t214 * t225 - t220 * t229) * t205 + (t238 * t212 + t214 * t216) * t206, t203 * r_i_i_C(1) - t204 * r_i_i_C(2); t204 * r_i_i_C(1) + t203 * r_i_i_C(2) + t214 * qJD(3) + t215 * t212 + (t231 * t212 - t234 * t214) * qJD(1), t230, t207, (t212 * t225 + t220 * t227) * t205 + (t212 * t216 - t238 * t214) * t206, -t201 * r_i_i_C(1) + t202 * r_i_i_C(2); 0, 0, 0, t205 * t217 + (-t220 * t206 - t236) * qJD(4), (t205 * t228 + t211 * t226) * r_i_i_C(2) + (qJD(4) * t205 * t211 - t213 * t226) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end