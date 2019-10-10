% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:41
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRPR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:04
	% EndTime: 2019-10-09 23:41:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:04
	% EndTime: 2019-10-09 23:41:04
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
	% StartTime: 2019-10-09 23:41:04
	% EndTime: 2019-10-09 23:41:04
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->6), mult. (20->10), div. (0->0), fcn. (12->2), ass. (0->5)
	t10 = -pkin(1) + r_i_i_C(2);
	t9 = r_i_i_C(3) + qJ(2);
	t8 = cos(qJ(1));
	t7 = sin(qJ(1));
	t1 = [t8 * qJD(2) + (t10 * t8 - t9 * t7) * qJD(1), qJD(1) * t8, 0, 0, 0, 0; t7 * qJD(2) + (t10 * t7 + t9 * t8) * qJD(1), qJD(1) * t7, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:04
	% EndTime: 2019-10-09 23:41:04
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (13->9), mult. (28->12), div. (0->0), fcn. (18->2), ass. (0->7)
	t12 = r_i_i_C(2) + qJ(2);
	t8 = sin(qJ(1));
	t11 = qJD(1) * t8;
	t10 = -pkin(1) - r_i_i_C(3) - qJ(3);
	t9 = cos(qJ(1));
	t7 = qJD(1) * t9;
	t1 = [t9 * qJD(2) - t8 * qJD(3) + (t10 * t9 - t12 * t8) * qJD(1), t7, -t11, 0, 0, 0; t8 * qJD(2) + t9 * qJD(3) + (t10 * t8 + t12 * t9) * qJD(1), t11, t7, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:04
	% EndTime: 2019-10-09 23:41:05
	% DurationCPUTime: 0.16s
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
	t25 = pkin(7) + r_i_i_C(3) - qJ(2);
	t22 = -r_i_i_C(1) * t18 - r_i_i_C(2) * t20 - pkin(1) - qJ(3);
	t1 = [t21 * qJD(2) - t29 * t19 + (t19 * t25 + t21 * t22) * qJD(1), t17, -t28, (t18 * t28 - t20 * t26) * r_i_i_C(2) + (-t18 * t26 - t20 * t28) * r_i_i_C(1), 0, 0; t19 * qJD(2) + t29 * t21 + (t19 * t22 - t21 * t25) * qJD(1), t28, t17, (-t17 * t18 - t20 * t27) * r_i_i_C(2) + (t17 * t20 - t18 * t27) * r_i_i_C(1), 0, 0; 0, 0, 0, -t23, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:05
	% EndTime: 2019-10-09 23:41:06
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (64->28), mult. (194->42), div. (0->0), fcn. (145->6), ass. (0->18)
	t167 = sin(qJ(4));
	t169 = cos(qJ(4));
	t165 = sin(pkin(9));
	t166 = cos(pkin(9));
	t175 = r_i_i_C(1) * t166 - r_i_i_C(2) * t165 + pkin(4);
	t180 = r_i_i_C(3) + qJ(5);
	t171 = -(t167 * t180 + t169 * t175) * qJD(4) + t169 * qJD(5);
	t182 = -qJD(3) + t171;
	t168 = sin(qJ(1));
	t179 = qJD(1) * t168;
	t170 = cos(qJ(1));
	t164 = qJD(1) * t170;
	t178 = qJD(4) * t167;
	t176 = qJD(4) * t180;
	t174 = r_i_i_C(1) * t165 + r_i_i_C(2) * t166 + pkin(7) - qJ(2);
	t173 = -qJD(4) * t175 + qJD(5);
	t172 = -t167 * t175 + t169 * t180 - pkin(1) - qJ(3);
	t1 = [t170 * qJD(2) + t182 * t168 + (t168 * t174 + t170 * t172) * qJD(1), t164, -t179, (t170 * t176 - t175 * t179) * t169 + (t170 * t173 - t179 * t180) * t167, t169 * t179 + t170 * t178, 0; t168 * qJD(2) - t182 * t170 + (t168 * t172 - t170 * t174) * qJD(1), t179, t164, (t164 * t175 + t168 * t176) * t169 + (t164 * t180 + t168 * t173) * t167, -t164 * t169 + t168 * t178, 0; 0, 0, 0, t171, qJD(4) * t169, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:05
	% EndTime: 2019-10-09 23:41:06
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (179->48), mult. (319->73), div. (0->0), fcn. (254->8), ass. (0->37)
	t203 = cos(pkin(9)) * pkin(5) + pkin(4);
	t212 = cos(qJ(4));
	t229 = t212 * qJD(5);
	t210 = sin(qJ(4));
	t237 = r_i_i_C(3) + pkin(8) + qJ(5);
	t239 = t237 * t210;
	t244 = (t203 * t212 + t239) * qJD(4) + qJD(3) - t229;
	t207 = pkin(9) + qJ(6);
	t204 = sin(t207);
	t205 = cos(t207);
	t218 = r_i_i_C(1) * t205 - r_i_i_C(2) * t204 + t203;
	t242 = t218 * t212 + t239;
	t213 = cos(qJ(1));
	t231 = qJD(6) * t210;
	t224 = qJD(1) + t231;
	t241 = t213 * t224;
	t211 = sin(qJ(1));
	t223 = qJD(1) * t210 + qJD(6);
	t233 = qJD(4) * t212;
	t238 = t211 * t233 + t223 * t213;
	t235 = qJD(1) * t211;
	t206 = qJD(1) * t213;
	t234 = qJD(4) * t210;
	t232 = qJD(4) * t213;
	t230 = qJD(6) * t212;
	t225 = t237 * t212;
	t221 = pkin(5) * sin(pkin(9)) + pkin(7) - qJ(2);
	t220 = r_i_i_C(1) * t204 + r_i_i_C(2) * t205;
	t219 = t224 * t211;
	t217 = -t203 * t210 - pkin(1) - qJ(3) + t225;
	t216 = t223 * t211 - t212 * t232;
	t214 = qJD(5) * t210 - t220 * t230 + (-t218 * t210 + t225) * qJD(4);
	t202 = t204 * t219 - t238 * t205;
	t201 = t238 * t204 + t205 * t219;
	t200 = t204 * t241 + t216 * t205;
	t199 = t216 * t204 - t205 * t241;
	t1 = [t202 * r_i_i_C(1) + t201 * r_i_i_C(2) + t213 * qJD(2) - t244 * t211 + (t221 * t211 + t217 * t213) * qJD(1), t206, -t235, t214 * t213 - t235 * t242, t210 * t232 + t212 * t235, t199 * r_i_i_C(1) + t200 * r_i_i_C(2); -t200 * r_i_i_C(1) + t199 * r_i_i_C(2) + t211 * qJD(2) + t244 * t213 + (t217 * t211 - t221 * t213) * qJD(1), t235, t206, t242 * t206 + t214 * t211, -t212 * t206 + t211 * t234, -t201 * r_i_i_C(1) + t202 * r_i_i_C(2); 0, 0, 0, -qJD(4) * t242 + t220 * t231 + t229, t233, (t204 * t230 + t205 * t234) * r_i_i_C(2) + (t204 * t234 - t205 * t230) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end