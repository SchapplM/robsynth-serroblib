% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:21
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPPR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:21:36
	% EndTime: 2019-10-10 09:21:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:21:36
	% EndTime: 2019-10-10 09:21:36
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
	% StartTime: 2019-10-10 09:21:36
	% EndTime: 2019-10-10 09:21:36
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(7) + r_i_i_C(3);
	t18 = sin(qJ(1));
	t27 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t26 = qJD(1) * t20;
	t25 = qJD(2) * t18;
	t24 = qJD(2) * t20;
	t17 = sin(qJ(2));
	t19 = cos(qJ(2));
	t23 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19;
	t22 = -r_i_i_C(1) * t19 + r_i_i_C(2) * t17 - pkin(1);
	t21 = t23 * qJD(2);
	t1 = [t23 * t25 + (-t28 * t18 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0, 0; -t20 * t21 + (t22 * t18 + t28 * t20) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0, 0; 0, -t21, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:21:37
	% EndTime: 2019-10-10 09:21:37
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (44->20), mult. (134->34), div. (0->0), fcn. (92->4), ass. (0->15)
	t134 = sin(qJ(2));
	t136 = cos(qJ(2));
	t146 = r_i_i_C(3) + qJ(3);
	t148 = pkin(2) + r_i_i_C(1);
	t149 = t148 * t134 - t146 * t136;
	t150 = t149 * qJD(2) - t134 * qJD(3);
	t147 = pkin(7) + r_i_i_C(2);
	t135 = sin(qJ(1));
	t145 = qJD(1) * t135;
	t137 = cos(qJ(1));
	t144 = qJD(1) * t137;
	t143 = qJD(2) * t137;
	t141 = -t146 * t134 - t148 * t136;
	t139 = -pkin(1) + t141;
	t1 = [t150 * t135 + (-t147 * t135 + t139 * t137) * qJD(1), (-t146 * t143 + t148 * t145) * t134 + (-t146 * t145 + (-t148 * qJD(2) + qJD(3)) * t137) * t136, -t134 * t145 + t136 * t143, 0, 0, 0; -t150 * t137 + (t139 * t135 + t147 * t137) * qJD(1), -t149 * t144 + (t141 * qJD(2) + qJD(3) * t136) * t135, t135 * qJD(2) * t136 + t134 * t144, 0, 0, 0; 0, -t150, qJD(2) * t134, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:21:36
	% EndTime: 2019-10-10 09:21:37
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (59->24), mult. (168->36), div. (0->0), fcn. (115->4), ass. (0->15)
	t19 = sin(qJ(2));
	t21 = cos(qJ(2));
	t28 = pkin(2) + pkin(3) - r_i_i_C(2);
	t33 = r_i_i_C(1) + qJ(3);
	t34 = t28 * t19 - t33 * t21;
	t35 = t34 * qJD(2) - t19 * qJD(3);
	t20 = sin(qJ(1));
	t32 = qJD(1) * t20;
	t22 = cos(qJ(1));
	t31 = qJD(1) * t22;
	t30 = qJD(2) * t22;
	t27 = pkin(7) - r_i_i_C(3) - qJ(4);
	t26 = -t33 * t19 - t28 * t21;
	t24 = -pkin(1) + t26;
	t1 = [-t22 * qJD(4) + t35 * t20 + (-t27 * t20 + t24 * t22) * qJD(1), (t28 * t32 - t33 * t30) * t19 + (-t33 * t32 + (-t28 * qJD(2) + qJD(3)) * t22) * t21, -t19 * t32 + t21 * t30, -t31, 0, 0; -t20 * qJD(4) - t35 * t22 + (t24 * t20 + t27 * t22) * qJD(1), -t34 * t31 + (t26 * qJD(2) + qJD(3) * t21) * t20, t20 * qJD(2) * t21 + t19 * t31, -t32, 0, 0; 0, -t35, qJD(2) * t19, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:21:37
	% EndTime: 2019-10-10 09:21:38
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (95->26), mult. (282->37), div. (0->0), fcn. (208->6), ass. (0->20)
	t162 = sin(pkin(9));
	t163 = cos(pkin(9));
	t164 = sin(qJ(2));
	t166 = cos(qJ(2));
	t175 = r_i_i_C(1) * t163 - r_i_i_C(2) * t162 + pkin(4) + qJ(3);
	t177 = pkin(2) + pkin(3) + r_i_i_C(3) + qJ(5);
	t172 = t177 * t164 - t175 * t166;
	t168 = -t172 * qJD(2) + t164 * qJD(3) + t166 * qJD(5);
	t184 = t168 - (t162 * r_i_i_C(1) + t163 * r_i_i_C(2) - pkin(7) + qJ(4)) * qJD(1);
	t165 = sin(qJ(1));
	t182 = qJD(1) * t165;
	t167 = cos(qJ(1));
	t181 = qJD(1) * t167;
	t180 = qJD(2) * t164;
	t179 = qJD(2) * t166;
	t178 = qJD(2) * t167;
	t173 = -t175 * t164 - t177 * t166;
	t170 = -qJD(4) + (-pkin(1) + t173) * qJD(1);
	t169 = t173 * qJD(2) + qJD(3) * t166 - qJD(5) * t164;
	t1 = [-t184 * t165 + t170 * t167, t169 * t167 + t172 * t182, -t164 * t182 + t166 * t178, -t181, -t164 * t178 - t166 * t182, 0; t170 * t165 + t184 * t167, t169 * t165 - t172 * t181, t164 * t181 + t165 * t179, -t182, -t165 * t180 + t166 * t181, 0; 0, t168, t180, 0, t179, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:21:38
	% EndTime: 2019-10-10 09:21:38
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (210->50), mult. (407->73), div. (0->0), fcn. (317->8), ass. (0->37)
	t209 = sin(qJ(2));
	t211 = cos(qJ(2));
	t230 = t211 * qJD(5);
	t237 = qJ(3) + cos(pkin(9)) * pkin(5) + pkin(4);
	t227 = pkin(2) + pkin(3) + r_i_i_C(3) + pkin(8) + qJ(5);
	t241 = t227 * t209;
	t246 = (-t237 * t211 + t241) * qJD(2) + (pkin(5) * sin(pkin(9)) - pkin(7) + qJ(4)) * qJD(1) - t209 * qJD(3) - t230;
	t206 = pkin(9) + qJ(6);
	t204 = sin(t206);
	t205 = cos(t206);
	t218 = r_i_i_C(1) * t205 - r_i_i_C(2) * t204 + t237;
	t244 = -t218 * t211 + t241;
	t212 = cos(qJ(1));
	t226 = qJD(6) * t209 + qJD(1);
	t243 = t212 * t226;
	t225 = qJD(1) * t209 + qJD(6);
	t210 = sin(qJ(1));
	t233 = qJD(2) * t211;
	t229 = t210 * t233;
	t239 = t225 * t212 + t229;
	t236 = qJD(1) * t210;
	t235 = qJD(1) * t212;
	t234 = qJD(2) * t209;
	t232 = qJD(2) * t212;
	t231 = qJD(6) * t211;
	t228 = t211 * t232;
	t222 = t227 * t211;
	t220 = t226 * t210;
	t217 = qJD(3) + (-r_i_i_C(1) * t204 - r_i_i_C(2) * t205) * qJD(6);
	t216 = t225 * t210 - t228;
	t215 = -qJD(4) + (-t237 * t209 - pkin(1) - t222) * qJD(1);
	t213 = -qJD(5) * t209 + t217 * t211 + (-t218 * t209 - t222) * qJD(2);
	t202 = t204 * t220 - t239 * t205;
	t201 = t239 * t204 + t205 * t220;
	t200 = t204 * t243 + t216 * t205;
	t199 = t216 * t204 - t205 * t243;
	t1 = [t202 * r_i_i_C(1) + t201 * r_i_i_C(2) + t246 * t210 + t215 * t212, t213 * t212 + t244 * t236, -t209 * t236 + t228, -t235, -t209 * t232 - t211 * t236, t199 * r_i_i_C(1) + t200 * r_i_i_C(2); -t200 * r_i_i_C(1) + t199 * r_i_i_C(2) + t215 * t210 - t246 * t212, t213 * t210 - t235 * t244, t209 * t235 + t229, -t236, -t210 * t234 + t211 * t235, -t201 * r_i_i_C(1) + t202 * r_i_i_C(2); 0, -qJD(2) * t244 + t217 * t209 + t230, t234, 0, t233, (-t204 * t231 - t205 * t234) * r_i_i_C(2) + (-t204 * t234 + t205 * t231) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end