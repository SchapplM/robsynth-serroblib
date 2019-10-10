% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:34
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRP4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:34:14
	% EndTime: 2019-10-10 00:34:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:34:14
	% EndTime: 2019-10-10 00:34:14
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
	% StartTime: 2019-10-10 00:34:14
	% EndTime: 2019-10-10 00:34:14
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
	% StartTime: 2019-10-10 00:34:14
	% EndTime: 2019-10-10 00:34:15
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (43->18), mult. (68->31), div. (0->0), fcn. (42->6), ass. (0->13)
	t24 = sin(qJ(3));
	t25 = cos(qJ(3));
	t26 = (r_i_i_C(1) * t24 + r_i_i_C(2) * t25) * qJD(3);
	t33 = pkin(7) + r_i_i_C(3);
	t32 = qJD(1) * t24;
	t31 = qJD(1) * t25;
	t30 = qJD(3) * t24;
	t29 = qJD(3) * t25;
	t27 = -r_i_i_C(1) * t25 + r_i_i_C(2) * t24 - pkin(2);
	t23 = qJ(1) + pkin(9);
	t22 = cos(t23);
	t21 = sin(t23);
	t1 = [t21 * t26 + (-cos(qJ(1)) * pkin(1) - t33 * t21 + t27 * t22) * qJD(1), 0, (t21 * t31 + t22 * t30) * r_i_i_C(2) + (t21 * t32 - t22 * t29) * r_i_i_C(1), 0, 0, 0; -t22 * t26 + (-sin(qJ(1)) * pkin(1) + t33 * t22 + t27 * t21) * qJD(1), 0, (t21 * t30 - t22 * t31) * r_i_i_C(2) + (-t21 * t29 - t22 * t32) * r_i_i_C(1), 0, 0, 0; 0, 0, -t26, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:34:15
	% EndTime: 2019-10-10 00:34:16
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (92->20), mult. (138->30), div. (0->0), fcn. (94->6), ass. (0->17)
	t147 = sin(qJ(3));
	t148 = cos(qJ(3));
	t158 = r_i_i_C(3) + qJ(4);
	t160 = pkin(3) - r_i_i_C(2);
	t153 = t160 * t147 - t158 * t148;
	t162 = qJD(1) * t153;
	t161 = t153 * qJD(3) - t147 * qJD(4);
	t159 = pkin(7) + r_i_i_C(1);
	t157 = qJD(1) * t147;
	t156 = qJD(3) * t148;
	t154 = -t158 * t147 - t160 * t148;
	t151 = -pkin(2) + t154;
	t150 = t154 * qJD(3) + qJD(4) * t148;
	t146 = qJ(1) + pkin(9);
	t145 = cos(t146);
	t144 = sin(t146);
	t1 = [t161 * t144 + (-cos(qJ(1)) * pkin(1) - t159 * t144 + t151 * t145) * qJD(1), 0, t144 * t162 + t150 * t145, -t144 * t157 + t145 * t156, 0, 0; -t161 * t145 + (-sin(qJ(1)) * pkin(1) + t159 * t145 + t151 * t144) * qJD(1), 0, t150 * t144 - t145 * t162, t144 * t156 + t145 * t157, 0, 0; 0, 0, -t161, qJD(3) * t147, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:34:16
	% EndTime: 2019-10-10 00:34:16
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (201->42), mult. (322->69), div. (0->0), fcn. (250->8), ass. (0->33)
	t204 = cos(qJ(3));
	t201 = sin(qJ(5));
	t203 = cos(qJ(5));
	t209 = r_i_i_C(1) * t201 + r_i_i_C(2) * t203 + qJ(4);
	t202 = sin(qJ(3));
	t219 = pkin(3) + pkin(8) + r_i_i_C(3);
	t229 = t219 * t202;
	t231 = -t209 * t204 + t229;
	t234 = qJD(1) * t231;
	t233 = (-qJ(4) * t204 + t229) * qJD(3) - t202 * qJD(4);
	t213 = qJD(5) * t202 + qJD(1);
	t222 = qJD(3) * t204;
	t228 = t213 * t201 - t203 * t222;
	t227 = t201 * t222 + t213 * t203;
	t226 = pkin(4) + pkin(7);
	t224 = qJD(1) * t202;
	t223 = qJD(3) * t202;
	t221 = qJD(5) * t204;
	t215 = t219 * t204;
	t212 = -qJD(5) - t224;
	t211 = t212 * t201;
	t210 = t212 * t203;
	t208 = -qJ(4) * t202 - pkin(2) - t215;
	t207 = qJD(4) + (r_i_i_C(1) * t203 - r_i_i_C(2) * t201) * qJD(5);
	t205 = t207 * t204 + (-t209 * t202 - t215) * qJD(3);
	t200 = qJ(1) + pkin(9);
	t199 = cos(t200);
	t198 = sin(t200);
	t197 = t198 * t211 + t227 * t199;
	t196 = t198 * t210 - t228 * t199;
	t195 = -t227 * t198 + t199 * t211;
	t194 = t228 * t198 + t199 * t210;
	t1 = [t195 * r_i_i_C(1) + t194 * r_i_i_C(2) + t233 * t198 + (-cos(qJ(1)) * pkin(1) - t226 * t198 + t208 * t199) * qJD(1), 0, t198 * t234 + t205 * t199, -t198 * t224 + t199 * t222, t196 * r_i_i_C(1) - t197 * r_i_i_C(2), 0; t197 * r_i_i_C(1) + t196 * r_i_i_C(2) - t233 * t199 + (-sin(qJ(1)) * pkin(1) + t226 * t199 + t208 * t198) * qJD(1), 0, t205 * t198 - t199 * t234, t198 * t222 + t199 * t224, -t194 * r_i_i_C(1) + t195 * r_i_i_C(2), 0; 0, 0, -qJD(3) * t231 + t207 * t202, t223, (-t201 * t223 + t203 * t221) * r_i_i_C(2) + (t201 * t221 + t203 * t223) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:34:16
	% EndTime: 2019-10-10 00:34:17
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (354->54), mult. (562->83), div. (0->0), fcn. (466->8), ass. (0->40)
	t241 = sin(qJ(3));
	t243 = cos(qJ(3));
	t240 = sin(qJ(5));
	t242 = cos(qJ(5));
	t270 = r_i_i_C(3) + qJ(6);
	t271 = r_i_i_C(1) + pkin(5);
	t273 = t271 * t240 - t270 * t242 + qJ(4);
	t246 = t273 * t243;
	t259 = pkin(3) + pkin(8) + r_i_i_C(2);
	t279 = qJD(1) * (-t259 * t241 + t246);
	t262 = qJD(3) * t243;
	t274 = t259 * qJD(3) + qJD(6) * t242 - qJD(4);
	t277 = -qJ(4) * t262 + t274 * t241;
	t239 = qJ(1) + pkin(9);
	t237 = sin(t239);
	t238 = cos(t239);
	t264 = qJD(1) * t241;
	t275 = t237 * t262 + t238 * t264;
	t272 = pkin(4) + pkin(7);
	t269 = t237 * t240;
	t268 = t237 * t242;
	t267 = t238 * t240;
	t266 = t238 * t242;
	t265 = t240 * t241;
	t263 = qJD(3) * t241;
	t261 = qJD(5) * t243;
	t260 = t240 * qJD(6);
	t255 = t238 * t262;
	t254 = qJD(5) * t241 + qJD(1);
	t253 = qJD(5) + t264;
	t251 = t238 * t265 + t268;
	t250 = -qJ(4) * t241 - t259 * t243 - pkin(2);
	t247 = t240 * t262 + t254 * t242;
	t245 = (t270 * t240 + t271 * t242) * qJD(5) - t274;
	t244 = t245 * t243 - t263 * t273;
	t232 = t247 * t238 - t253 * t269;
	t231 = -t242 * t255 + t251 * qJD(5) + (t241 * t268 + t267) * qJD(1);
	t230 = t247 * t237 + t253 * t267;
	t229 = -qJD(5) * t266 - t275 * t242 + t254 * t269;
	t1 = [t238 * t260 - t271 * t230 - t270 * t229 + t277 * t237 + (-cos(qJ(1)) * pkin(1) - t272 * t237 + t250 * t238) * qJD(1), 0, -t237 * t279 + t244 * t238, -t237 * t264 + t255, t251 * qJD(6) - t271 * t231 + t270 * t232, t231; t237 * t260 + t271 * t232 + t270 * t231 - t277 * t238 + (-sin(qJ(1)) * pkin(1) + t272 * t238 + t250 * t237) * qJD(1), 0, t244 * t237 + t238 * t279, t275, -(-t237 * t265 + t266) * qJD(6) + t270 * t230 - t271 * t229, t229; 0, 0, qJD(3) * t246 + t245 * t241, t263, (-t270 * t261 + t271 * t263) * t242 + (t270 * t263 + (t271 * qJD(5) - qJD(6)) * t243) * t240, -t240 * t261 - t242 * t263;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end