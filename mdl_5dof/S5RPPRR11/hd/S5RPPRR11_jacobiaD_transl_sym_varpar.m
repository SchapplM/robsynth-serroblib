% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPRR11
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPPRR11_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPRR11_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR11_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_jacobiaD_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:06:08
	% EndTime: 2019-12-31 18:06:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:06:08
	% EndTime: 2019-12-31 18:06:08
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
	% StartTime: 2019-12-31 18:06:08
	% EndTime: 2019-12-31 18:06:08
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
	% StartTime: 2019-12-31 18:06:08
	% EndTime: 2019-12-31 18:06:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (13->9), mult. (28->12), div. (0->0), fcn. (18->2), ass. (0->7)
	t12 = r_i_i_C(2) + qJ(2);
	t8 = sin(qJ(1));
	t11 = qJD(1) * t8;
	t10 = -pkin(1) - r_i_i_C(3) - qJ(3);
	t9 = cos(qJ(1));
	t7 = qJD(1) * t9;
	t1 = [t9 * qJD(2) - t8 * qJD(3) + (t10 * t9 - t12 * t8) * qJD(1), t7, -t11, 0, 0; t8 * qJD(2) + t9 * qJD(3) + (t10 * t8 + t12 * t9) * qJD(1), t11, t7, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:06:08
	% EndTime: 2019-12-31 18:06:08
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (28->20), mult. (80->31), div. (0->0), fcn. (52->4), ass. (0->13)
	t18 = sin(qJ(4));
	t20 = cos(qJ(4));
	t23 = (r_i_i_C(1) * t20 - r_i_i_C(2) * t18) * qJD(4);
	t29 = qJD(3) + t23;
	t19 = sin(qJ(1));
	t28 = qJD(1) * t19;
	t21 = cos(qJ(1));
	t17 = qJD(1) * t21;
	t27 = qJD(4) * t19;
	t26 = qJD(4) * t21;
	t25 = pkin(6) + r_i_i_C(3) - qJ(2);
	t22 = -r_i_i_C(1) * t18 - r_i_i_C(2) * t20 - pkin(1) - qJ(3);
	t1 = [t21 * qJD(2) - t29 * t19 + (t25 * t19 + t22 * t21) * qJD(1), t17, -t28, (t18 * t28 - t20 * t26) * r_i_i_C(2) + (-t18 * t26 - t20 * t28) * r_i_i_C(1), 0; t19 * qJD(2) + t29 * t21 + (t22 * t19 - t25 * t21) * qJD(1), t28, t17, (-t18 * t17 - t20 * t27) * r_i_i_C(2) + (t20 * t17 - t18 * t27) * r_i_i_C(1), 0; 0, 0, 0, -t23, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:06:09
	% EndTime: 2019-12-31 18:06:09
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (92->44), mult. (286->73), div. (0->0), fcn. (223->6), ass. (0->33)
	t206 = cos(qJ(4));
	t203 = sin(qJ(4));
	t228 = pkin(7) + r_i_i_C(3);
	t229 = t228 * t203;
	t233 = (pkin(4) * t206 + t229) * qJD(4) + qJD(3);
	t202 = sin(qJ(5));
	t205 = cos(qJ(5));
	t212 = r_i_i_C(1) * t205 - r_i_i_C(2) * t202 + pkin(4);
	t231 = t212 * t206 + t229;
	t226 = pkin(6) - qJ(2);
	t207 = cos(qJ(1));
	t225 = t205 * t207;
	t204 = sin(qJ(1));
	t224 = qJD(1) * t204;
	t201 = qJD(1) * t207;
	t223 = qJD(4) * t203;
	t222 = qJD(4) * t206;
	t221 = qJD(4) * t207;
	t220 = qJD(5) * t203;
	t219 = qJD(5) * t206;
	t216 = t228 * t206;
	t215 = qJD(1) + t220;
	t214 = qJD(1) * t203 + qJD(5);
	t213 = r_i_i_C(1) * t202 + r_i_i_C(2) * t205;
	t211 = t215 * t202;
	t210 = t213 * qJD(5);
	t209 = -pkin(4) * t203 - pkin(1) - qJ(3) + t216;
	t208 = t214 * t204 - t206 * t221;
	t200 = -t214 * t225 + (-t205 * t222 + t211) * t204;
	t199 = t215 * t205 * t204 + (t204 * t222 + t214 * t207) * t202;
	t198 = t208 * t205 + t207 * t211;
	t197 = t208 * t202 - t215 * t225;
	t1 = [t200 * r_i_i_C(1) + t199 * r_i_i_C(2) + t207 * qJD(2) - t233 * t204 + (t226 * t204 + t209 * t207) * qJD(1), t201, -t224, (-t212 * t221 - t228 * t224) * t203 + (-t212 * t224 + (t228 * qJD(4) - t210) * t207) * t206, t197 * r_i_i_C(1) + t198 * r_i_i_C(2); -t198 * r_i_i_C(1) + t197 * r_i_i_C(2) + t204 * qJD(2) + t233 * t207 + (t209 * t204 - t226 * t207) * qJD(1), t224, t201, t231 * t201 + (-t206 * t210 + (-t212 * t203 + t216) * qJD(4)) * t204, -t199 * r_i_i_C(1) + t200 * r_i_i_C(2); 0, 0, 0, -t231 * qJD(4) + t213 * t220, (t202 * t219 + t205 * t223) * r_i_i_C(2) + (t202 * t223 - t205 * t219) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end