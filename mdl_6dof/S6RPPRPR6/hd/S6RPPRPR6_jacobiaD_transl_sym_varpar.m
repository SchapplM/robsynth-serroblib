% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRPR6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:42
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRPR6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:44
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
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:44
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
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:44
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
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:45
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
	% StartTime: 2019-10-09 23:42:45
	% EndTime: 2019-10-09 23:42:45
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (53->25), mult. (150->38), div. (0->0), fcn. (104->4), ass. (0->16)
	t139 = sin(qJ(4));
	t141 = cos(qJ(4));
	t151 = r_i_i_C(3) + qJ(5);
	t152 = pkin(4) - r_i_i_C(2);
	t143 = -(t151 * t139 + t152 * t141) * qJD(4) + t141 * qJD(5);
	t154 = -qJD(3) + t143;
	t140 = sin(qJ(1));
	t150 = qJD(1) * t140;
	t142 = cos(qJ(1));
	t138 = qJD(1) * t142;
	t149 = qJD(4) * t139;
	t147 = pkin(7) + r_i_i_C(1) - qJ(2);
	t146 = qJD(4) * t151;
	t145 = -t152 * qJD(4) + qJD(5);
	t144 = -t152 * t139 + t151 * t141 - pkin(1) - qJ(3);
	t1 = [t142 * qJD(2) + t154 * t140 + (t147 * t140 + t144 * t142) * qJD(1), t138, -t150, (t142 * t146 - t152 * t150) * t141 + (t145 * t142 - t151 * t150) * t139, t141 * t150 + t142 * t149, 0; t140 * qJD(2) - t154 * t142 + (t144 * t140 - t147 * t142) * qJD(1), t150, t138, (t152 * t138 + t140 * t146) * t141 + (t151 * t138 + t145 * t140) * t139, -t141 * t138 + t140 * t149, 0; 0, 0, 0, t143, qJD(4) * t141, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:46
	% EndTime: 2019-10-09 23:42:46
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (110->48), mult. (334->77), div. (0->0), fcn. (260->6), ass. (0->33)
	t198 = sin(qJ(4));
	t201 = cos(qJ(4));
	t216 = pkin(4) + pkin(8) + r_i_i_C(3);
	t212 = t216 * t201;
	t226 = (qJ(5) * t198 + t212) * qJD(4) - t201 * qJD(5) + qJD(3);
	t202 = cos(qJ(1));
	t220 = qJD(1) * t201;
	t210 = qJD(6) + t220;
	t224 = t210 * t202;
	t199 = sin(qJ(1));
	t219 = qJD(4) * t198;
	t213 = t202 * t219;
	t223 = t199 * t210 + t213;
	t221 = qJD(1) * t199;
	t196 = qJD(1) * t202;
	t218 = qJD(4) * t201;
	t217 = qJD(6) * t198;
	t215 = pkin(5) + pkin(7) - qJ(2);
	t214 = t199 * t219;
	t211 = qJD(6) * t201 + qJD(1);
	t208 = t211 * t202;
	t197 = sin(qJ(6));
	t200 = cos(qJ(6));
	t207 = r_i_i_C(1) * t197 + r_i_i_C(2) * t200 + qJ(5);
	t206 = qJD(4) * t207;
	t205 = qJD(5) + (r_i_i_C(1) * t200 - r_i_i_C(2) * t197) * qJD(6);
	t204 = qJ(5) * t201 - t198 * t216 - pkin(1) - qJ(3);
	t203 = -qJD(4) * t216 + t205;
	t195 = -t197 * t223 + t200 * t208;
	t194 = t197 * t208 + t200 * t223;
	t193 = t211 * t200 * t199 + (-t214 + t224) * t197;
	t192 = -t200 * t224 + (t197 * t211 + t200 * t219) * t199;
	t1 = [t193 * r_i_i_C(1) - t192 * r_i_i_C(2) + t202 * qJD(2) - t226 * t199 + (t199 * t215 + t202 * t204) * qJD(1), t196, -t221, (t202 * t206 - t216 * t221) * t201 + (t202 * t203 - t207 * t221) * t198, t199 * t220 + t213, r_i_i_C(1) * t194 + r_i_i_C(2) * t195; -t195 * r_i_i_C(1) + t194 * r_i_i_C(2) + t199 * qJD(2) + t226 * t202 + (t199 * t204 - t202 * t215) * qJD(1), t221, t196, (t196 * t216 + t199 * t206) * t201 + (t196 * t207 + t199 * t203) * t198, -t196 * t201 + t214, r_i_i_C(1) * t192 + r_i_i_C(2) * t193; 0, 0, 0, t205 * t201 + (-t198 * t207 - t212) * qJD(4), t218, (-t197 * t218 - t200 * t217) * r_i_i_C(2) + (-t197 * t217 + t200 * t218) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end