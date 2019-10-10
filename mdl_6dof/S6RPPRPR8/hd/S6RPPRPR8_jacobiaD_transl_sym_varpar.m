% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:46
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRPR8_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR8_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR8_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:46:04
	% EndTime: 2019-10-09 23:46:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:46:04
	% EndTime: 2019-10-09 23:46:04
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
	% StartTime: 2019-10-09 23:46:04
	% EndTime: 2019-10-09 23:46:05
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
	% StartTime: 2019-10-09 23:46:04
	% EndTime: 2019-10-09 23:46:04
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (15->10), mult. (36->14), div. (0->0), fcn. (24->4), ass. (0->7)
	t14 = sin(qJ(1));
	t18 = qJD(1) * t14;
	t17 = -pkin(1) - r_i_i_C(3) - qJ(3);
	t16 = r_i_i_C(1) * sin(pkin(9)) + r_i_i_C(2) * cos(pkin(9)) + qJ(2);
	t15 = cos(qJ(1));
	t11 = qJD(1) * t15;
	t1 = [t15 * qJD(2) - t14 * qJD(3) + (-t16 * t14 + t17 * t15) * qJD(1), t11, -t18, 0, 0, 0; t14 * qJD(2) + t15 * qJD(3) + (t17 * t14 + t16 * t15) * qJD(1), t18, t11, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:46:04
	% EndTime: 2019-10-09 23:46:05
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (48->21), mult. (82->32), div. (0->0), fcn. (54->5), ass. (0->14)
	t24 = pkin(9) + qJ(4);
	t21 = sin(t24);
	t22 = cos(t24);
	t36 = (r_i_i_C(1) * t22 - r_i_i_C(2) * t21) * qJD(4);
	t27 = sin(qJ(1));
	t35 = qJD(1) * t27;
	t28 = cos(qJ(1));
	t23 = qJD(1) * t28;
	t34 = qJD(4) * t27;
	t33 = qJD(4) * t28;
	t32 = -pkin(1) - r_i_i_C(3) - pkin(7) - qJ(3);
	t30 = pkin(3) * sin(pkin(9)) + r_i_i_C(1) * t21 + r_i_i_C(2) * t22 + qJ(2);
	t29 = qJD(2) + t36;
	t1 = [-t27 * qJD(3) + t29 * t28 + (-t30 * t27 + t32 * t28) * qJD(1), t23, -t35, (-t21 * t23 - t22 * t34) * r_i_i_C(2) + (-t21 * t34 + t22 * t23) * r_i_i_C(1), 0, 0; t28 * qJD(3) + t29 * t27 + (t32 * t27 + t30 * t28) * qJD(1), t35, t23, (-t21 * t35 + t22 * t33) * r_i_i_C(2) + (t21 * t33 + t22 * t35) * r_i_i_C(1), 0, 0; 0, 0, 0, -t36, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:46:05
	% EndTime: 2019-10-09 23:46:05
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (101->26), mult. (152->41), div. (0->0), fcn. (106->5), ass. (0->19)
	t145 = pkin(9) + qJ(4);
	t142 = sin(t145);
	t143 = cos(t145);
	t159 = r_i_i_C(3) + qJ(5);
	t160 = pkin(4) - r_i_i_C(2);
	t164 = -(t159 * t142 + t160 * t143) * qJD(4) + t143 * qJD(5);
	t163 = t160 * qJD(4) - qJD(5);
	t161 = t160 * t142 - t159 * t143 + pkin(3) * sin(pkin(9)) + qJ(2);
	t148 = sin(qJ(1));
	t158 = qJD(1) * t148;
	t149 = cos(qJ(1));
	t144 = qJD(1) * t149;
	t157 = qJD(4) * t148;
	t156 = qJD(4) * t149;
	t154 = -pkin(1) - r_i_i_C(1) - pkin(7) - qJ(3);
	t153 = qJD(1) * t160;
	t151 = qJD(1) * t159;
	t150 = qJD(2) - t164;
	t1 = [-t148 * qJD(3) + t150 * t149 + (-t161 * t148 + t154 * t149) * qJD(1), t144, -t158, (t149 * t153 + t159 * t157) * t143 + (-t163 * t148 + t149 * t151) * t142, t142 * t157 - t143 * t144, 0; qJD(3) * t149 + t150 * t148 + (t154 * t148 + t161 * t149) * qJD(1), t158, t144, (t148 * t153 - t159 * t156) * t143 + (t148 * t151 + t163 * t149) * t142, -t142 * t156 - t143 * t158, 0; 0, 0, 0, t164, qJD(4) * t143, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:46:06
	% EndTime: 2019-10-09 23:46:06
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (202->49), mult. (336->81), div. (0->0), fcn. (262->7), ass. (0->36)
	t208 = pkin(9) + qJ(4);
	t206 = cos(t208);
	t229 = pkin(4) + pkin(8) + r_i_i_C(3);
	t239 = t229 * t206;
	t205 = sin(t208);
	t238 = t229 * t205 + pkin(3) * sin(pkin(9)) - qJ(5) * t206 + qJ(2);
	t211 = sin(qJ(6));
	t213 = cos(qJ(6));
	t217 = qJD(5) + (r_i_i_C(1) * t213 - r_i_i_C(2) * t211) * qJD(6);
	t237 = -t229 * qJD(4) + t217;
	t214 = cos(qJ(1));
	t236 = t213 * t214;
	t212 = sin(qJ(1));
	t235 = qJD(1) * t212;
	t207 = qJD(1) * t214;
	t234 = qJD(4) * t206;
	t233 = qJD(4) * t212;
	t232 = qJD(4) * t213;
	t231 = qJD(4) * t214;
	t230 = qJD(6) * t205;
	t228 = -pkin(1) - pkin(5) - pkin(7) - qJ(3);
	t227 = t205 * t233;
	t226 = t205 * t231;
	t225 = qJD(1) * t229;
	t224 = qJD(6) * t206 + qJD(1);
	t223 = qJD(1) * t206 + qJD(6);
	t221 = t224 * t211;
	t220 = r_i_i_C(1) * t211 + r_i_i_C(2) * t213 + qJ(5);
	t218 = qJD(1) * t220;
	t216 = t223 * t212 + t226;
	t215 = -t206 * qJD(5) + qJD(2) + (qJ(5) * t205 + t239) * qJD(4);
	t204 = t216 * t211 - t224 * t236;
	t203 = t216 * t213 + t214 * t221;
	t202 = t224 * t213 * t212 + (t223 * t214 - t227) * t211;
	t201 = -t223 * t236 + (t205 * t232 + t221) * t212;
	t1 = [t204 * r_i_i_C(1) + t203 * r_i_i_C(2) - t212 * qJD(3) + t215 * t214 + (-t238 * t212 + t228 * t214) * qJD(1), t207, -t235, (t214 * t225 + t220 * t233) * t206 + (t237 * t212 + t214 * t218) * t205, -t206 * t207 + t227, t201 * r_i_i_C(1) + t202 * r_i_i_C(2); -t202 * r_i_i_C(1) + t201 * r_i_i_C(2) + t214 * qJD(3) + t215 * t212 + (t228 * t212 + t238 * t214) * qJD(1), t235, t207, (t212 * t225 - t220 * t231) * t206 + (t212 * t218 - t237 * t214) * t205, -t206 * t235 - t226, -t203 * r_i_i_C(1) + t204 * r_i_i_C(2); 0, 0, 0, t217 * t206 + (-t220 * t205 - t239) * qJD(4), t234, (-t211 * t234 - t213 * t230) * r_i_i_C(2) + (t206 * t232 - t211 * t230) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end