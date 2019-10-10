% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPPR4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:20
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
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
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (10->7), mult. (28->12), div. (0->0), fcn. (18->4), ass. (0->5)
	t16 = r_i_i_C(3) + qJ(2);
	t15 = -r_i_i_C(1) * cos(pkin(9)) + r_i_i_C(2) * sin(pkin(9)) - pkin(1);
	t14 = cos(qJ(1));
	t13 = sin(qJ(1));
	t1 = [t14 * qJD(2) + (-t16 * t13 + t15 * t14) * qJD(1), qJD(1) * t14, 0, 0, 0, 0; t13 * qJD(2) + (t15 * t13 + t16 * t14) * qJD(1), qJD(1) * t13, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (43->19), mult. (70->32), div. (0->0), fcn. (46->5), ass. (0->14)
	t35 = r_i_i_C(3) + pkin(7) + qJ(2);
	t26 = sin(qJ(1));
	t34 = qJD(1) * t26;
	t27 = cos(qJ(1));
	t33 = qJD(1) * t27;
	t32 = qJD(3) * t26;
	t31 = qJD(3) * t27;
	t24 = pkin(9) + qJ(3);
	t22 = sin(t24);
	t23 = cos(t24);
	t30 = r_i_i_C(1) * t22 + r_i_i_C(2) * t23;
	t29 = -r_i_i_C(1) * t23 + r_i_i_C(2) * t22 - cos(pkin(9)) * pkin(2) - pkin(1);
	t28 = t30 * qJD(3);
	t1 = [t27 * qJD(2) + t30 * t32 + (-t35 * t26 + t29 * t27) * qJD(1), t33, (t22 * t31 + t23 * t34) * r_i_i_C(2) + (t22 * t34 - t23 * t31) * r_i_i_C(1), 0, 0, 0; t26 * qJD(2) - t27 * t28 + (t29 * t26 + t35 * t27) * qJD(1), t34, (t22 * t32 - t23 * t33) * r_i_i_C(2) + (-t22 * t33 - t23 * t32) * r_i_i_C(1), 0, 0, 0; 0, 0, -t28, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:29
	% EndTime: 2019-10-10 00:20:30
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (116->28), mult. (184->43), div. (0->0), fcn. (139->7), ass. (0->18)
	t172 = pkin(9) + qJ(3);
	t170 = sin(t172);
	t171 = cos(t172);
	t173 = sin(pkin(10));
	t174 = cos(pkin(10));
	t182 = r_i_i_C(1) * t174 - r_i_i_C(2) * t173 + pkin(3);
	t188 = r_i_i_C(3) + qJ(4);
	t190 = (t182 * t170 - t188 * t171) * qJD(3) - t170 * qJD(4);
	t176 = sin(qJ(1));
	t187 = qJD(1) * t176;
	t177 = cos(qJ(1));
	t186 = qJD(1) * t177;
	t185 = qJD(3) * t171;
	t183 = qJD(3) * t188;
	t181 = t173 * r_i_i_C(1) + t174 * r_i_i_C(2) + pkin(7) + qJ(2);
	t180 = -t182 * qJD(3) + qJD(4);
	t179 = -t188 * t170 - t182 * t171 - cos(pkin(9)) * pkin(2) - pkin(1);
	t1 = [t177 * qJD(2) + t190 * t176 + (-t181 * t176 + t179 * t177) * qJD(1), t186, (-t177 * t183 + t182 * t187) * t170 + (t180 * t177 - t188 * t187) * t171, -t170 * t187 + t177 * t185, 0, 0; t176 * qJD(2) - t190 * t177 + (t179 * t176 + t181 * t177) * qJD(1), t187, (-t176 * t183 - t182 * t186) * t170 + (t180 * t176 + t188 * t186) * t171, t170 * t186 + t176 * t185, 0, 0; 0, 0, -t190, qJD(3) * t170, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:29
	% EndTime: 2019-10-10 00:20:30
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (177->42), mult. (292->65), div. (0->0), fcn. (235->7), ass. (0->34)
	t189 = pkin(9) + qJ(3);
	t187 = sin(t189);
	t190 = sin(pkin(10));
	t191 = cos(pkin(10));
	t215 = r_i_i_C(3) + qJ(5);
	t218 = pkin(4) + r_i_i_C(1);
	t219 = t215 * t190 + t218 * t191 + pkin(3);
	t188 = cos(t189);
	t216 = r_i_i_C(2) + qJ(4);
	t220 = t216 * t188;
	t223 = t219 * t187 - t220;
	t206 = qJD(5) * t190;
	t199 = t187 * qJD(4) + t188 * t206;
	t222 = (-pkin(3) * t187 + t220) * qJD(3) + t199;
	t193 = sin(qJ(1));
	t214 = t190 * t193;
	t194 = cos(qJ(1));
	t213 = t190 * t194;
	t212 = t193 * t191;
	t211 = t194 * t191;
	t210 = qJD(1) * t193;
	t209 = qJD(1) * t194;
	t208 = qJD(3) * t193;
	t207 = qJD(3) * t194;
	t205 = t187 * t208;
	t204 = t187 * t207;
	t203 = t216 * t187;
	t200 = -t191 * qJD(5) + qJD(2);
	t198 = -pkin(3) * t188 - cos(pkin(9)) * pkin(2) - pkin(1) - t203;
	t195 = -t187 * t206 + qJD(4) * t188 + (-t188 * t219 - t203) * qJD(3);
	t192 = -pkin(7) - qJ(2);
	t184 = -t190 * t205 + (t188 * t213 - t212) * qJD(1);
	t182 = t190 * t204 + (t188 * t214 + t211) * qJD(1);
	t1 = [t200 * t194 + t218 * (t191 * t205 + (-t188 * t211 - t214) * qJD(1)) - t215 * t184 - t222 * t193 + (t193 * t192 + t198 * t194) * qJD(1), t209, t195 * t194 + t223 * t210, -t187 * t210 + t188 * t207, -t182, 0; t200 * t193 + t218 * (-t191 * t204 + (-t188 * t212 + t213) * qJD(1)) - t215 * t182 + t222 * t194 + (-t194 * t192 + t198 * t193) * qJD(1), t210, t195 * t193 - t209 * t223, t187 * t209 + t188 * t208, t184, 0; 0, 0, -qJD(3) * t223 + t199, qJD(3) * t187, qJD(3) * t188 * t190, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:29
	% EndTime: 2019-10-10 00:20:30
	% DurationCPUTime: 0.64s
	% Computational Cost: add. (380->71), mult. (665->113), div. (0->0), fcn. (606->9), ass. (0->51)
	t258 = pkin(9) + qJ(3);
	t256 = sin(t258);
	t259 = sin(pkin(10));
	t260 = cos(pkin(10));
	t262 = sin(qJ(6));
	t264 = cos(qJ(6));
	t296 = pkin(4) + pkin(5);
	t271 = t264 * r_i_i_C(1) - t262 * r_i_i_C(2) + t296;
	t272 = t262 * r_i_i_C(1) + t264 * r_i_i_C(2) + qJ(5);
	t297 = t272 * t259 + t271 * t260 + pkin(3);
	t257 = cos(t258);
	t283 = -r_i_i_C(3) - pkin(8) + qJ(4);
	t299 = t283 * t257;
	t302 = t297 * t256 - t299;
	t285 = t256 * qJD(4);
	t301 = (-t256 * pkin(3) + t299) * qJD(3) + t285;
	t273 = t259 * t262 + t260 * t264;
	t274 = t259 * t264 - t260 * t262;
	t269 = t274 * r_i_i_C(1) - t273 * r_i_i_C(2);
	t298 = t259 * qJD(5) + t269 * qJD(6);
	t265 = cos(qJ(1));
	t294 = t260 * t265;
	t263 = sin(qJ(1));
	t293 = t263 * t259;
	t292 = t263 * t260;
	t291 = t265 * t259;
	t290 = qJD(1) * t263;
	t289 = qJD(1) * t265;
	t288 = qJD(3) * t257;
	t287 = qJD(3) * t263;
	t286 = qJD(3) * t265;
	t282 = t257 * t291;
	t281 = t256 * t287;
	t280 = t256 * t286;
	t279 = t283 * t256;
	t247 = t257 * t293 + t294;
	t248 = t257 * t292 - t291;
	t276 = -t247 * t264 + t248 * t262;
	t275 = t247 * t262 + t248 * t264;
	t250 = t257 * t294 + t293;
	t270 = -pkin(3) * t257 - cos(pkin(9)) * pkin(2) - pkin(1) - t279;
	t266 = t257 * qJD(4) - t298 * t256 + (-t257 * t297 - t279) * qJD(3);
	t261 = -pkin(7) - qJ(2);
	t249 = t282 - t292;
	t246 = t250 * qJD(1) - t260 * t281;
	t245 = qJD(1) * t282 - t259 * t281 - t260 * t290;
	t244 = -t248 * qJD(1) - t260 * t280;
	t243 = t247 * qJD(1) + t259 * t280;
	t242 = -t243 * t262 + t244 * t264 + (t249 * t264 - t250 * t262) * qJD(6);
	t241 = -t243 * t264 - t244 * t262 + (-t249 * t262 - t250 * t264) * qJD(6);
	t1 = [t265 * qJD(2) - t247 * qJD(5) - t272 * t245 - t271 * t246 + (t276 * r_i_i_C(1) + t275 * r_i_i_C(2)) * qJD(6) - t301 * t263 + (t263 * t261 + t270 * t265) * qJD(1), t289, t266 * t265 + t302 * t290, -t256 * t290 + t257 * t286, -t243, r_i_i_C(1) * t241 - r_i_i_C(2) * t242; t242 * r_i_i_C(1) + t241 * r_i_i_C(2) - t243 * qJ(5) + t263 * qJD(2) + t249 * qJD(5) + t296 * t244 + t301 * t265 + (-t261 * t265 + t270 * t263) * qJD(1), t290, t266 * t263 - t289 * t302, t256 * t289 + t257 * t287, t245, (t245 * t264 - t246 * t262) * r_i_i_C(1) + (-t245 * t262 - t246 * t264) * r_i_i_C(2) + (-t275 * r_i_i_C(1) + t276 * r_i_i_C(2)) * qJD(6); 0, 0, -qJD(3) * t302 + t298 * t257 + t285, qJD(3) * t256, t259 * t288, (-t273 * r_i_i_C(1) - t274 * r_i_i_C(2)) * t256 * qJD(6) + t269 * t288;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end