% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5PRRRR12_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PRRRR12_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR12_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
JaD_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (46->7), mult. (159->21), div. (18->4), fcn. (175->5), ass. (0->15)
	t39 = cos(pkin(11));
	t41 = sin(qJ(2));
	t42 = cos(qJ(2));
	t45 = sin(pkin(11)) * cos(pkin(5));
	t37 = t39 * t42 - t41 * t45;
	t34 = 0.1e1 / t37 ^ 2;
	t49 = qJD(2) * t34;
	t36 = t39 * t41 + t42 * t45;
	t33 = t36 ^ 2;
	t30 = t33 * t34 + 0.1e1;
	t46 = t37 * t49;
	t47 = t36 / t37 * t49;
	t48 = (t33 * t47 + t36 * t46) / t30 ^ 2;
	t28 = 0.1e1 / t30;
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, -0.2e1 * t48 + 0.2e1 * (t28 * t46 + (t28 * t47 - t34 * t48) * t36) * t36, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (728->15), mult. (516->22), div. (50->4), fcn. (460->8), ass. (0->17)
	t89 = qJ(2) + qJ(3);
	t84 = pkin(5) + t89;
	t85 = pkin(5) - t89;
	t90 = sin(pkin(11));
	t91 = cos(pkin(11));
	t74 = t91 * sin(t89) + t90 * (cos(t85) / 0.2e1 + cos(t84) / 0.2e1);
	t68 = t74 * (qJD(2) + qJD(3));
	t94 = t90 * (sin(t85) / 0.2e1 - sin(t84) / 0.2e1) + t91 * cos(t89);
	t72 = 0.1e1 / t94 ^ 2;
	t102 = t72 * t74 ^ 2;
	t71 = 0.1e1 / t94;
	t101 = t71 * t102;
	t100 = t72 * t94;
	t97 = t100 * t68;
	t67 = 0.1e1 + t102;
	t63 = 0.2e1 * (-t71 * t94 - t102) / t67 ^ 2 * (t68 * t101 + t97) + (0.2e1 * t97 - (-t100 + t71 - 0.2e1 * t101) * t68) / t67;
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, t63, t63, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (1587->17), mult. (774->27), div. (75->4), fcn. (690->8), ass. (0->22)
	t103 = cos(pkin(11));
	t102 = sin(pkin(11));
	t101 = qJ(2) + qJ(3) + qJ(4);
	t96 = pkin(5) + t101;
	t97 = pkin(5) - t101;
	t111 = t102 * (sin(t97) / 0.2e1 - sin(t96) / 0.2e1);
	t99 = cos(t101);
	t106 = t103 * t99 + t111;
	t84 = 0.1e1 / t106 ^ 2;
	t112 = t102 * (cos(t97) / 0.2e1 + cos(t96) / 0.2e1);
	t98 = sin(t101);
	t86 = t103 * t98 + t112;
	t114 = t84 * t86 ^ 2;
	t83 = 0.1e1 / t106;
	t113 = t83 * t114;
	t100 = qJD(2) + qJD(3) + qJD(4);
	t110 = t100 * t103;
	t109 = t86 * t84 * (t100 * t111 + t99 * t110);
	t80 = -t100 * t112 - t98 * t110;
	t79 = 0.1e1 + t114;
	t75 = 0.2e1 * (-t106 * t83 - t114) / t79 ^ 2 * (-t80 * t113 + t109) + (0.2e1 * t109 + (-t106 * t84 - 0.2e1 * t113 + t83) * t80) / t79;
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, t75, t75, t75, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:13
	% EndTime: 2024-09-28 18:09:13
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (10430->74), mult. (12726->168), div. (796->12), fcn. (16164->13), ass. (0->90)
	t281 = cos(pkin(6));
	t283 = sin(qJ(5));
	t277 = sin(pkin(11));
	t282 = cos(pkin(5));
	t315 = t277 * t282;
	t305 = t283 * t315;
	t280 = cos(pkin(11));
	t284 = cos(qJ(5));
	t311 = t280 * t284;
	t263 = t281 * t305 - t311;
	t304 = t284 * t315;
	t312 = t280 * t283;
	t265 = t281 * t312 + t304;
	t276 = qJ(2) + qJ(3) + qJ(4);
	t273 = sin(t276);
	t274 = cos(t276);
	t278 = sin(pkin(6));
	t279 = sin(pkin(5));
	t314 = t278 * t279;
	t306 = t277 * t314;
	t293 = -t263 * t274 - t265 * t273 + t283 * t306;
	t230 = 0.1e1 / t293 ^ 2;
	t292 = t281 * t304 + t312;
	t325 = t281 * t311 - t305;
	t235 = -t273 * t325 - t274 * t292 + t284 * t306;
	t232 = t235 ^ 2;
	t327 = t232 * t230;
	t275 = qJD(2) + qJD(3) + qJD(4);
	t326 = t265 * qJD(5) + t275 * t292;
	t294 = t274 * t278 * t282 + t279 * t281;
	t318 = t273 * t277;
	t251 = -t278 * t318 + t294 * t280;
	t303 = t274 * t314;
	t262 = t282 * t281 - t303;
	t242 = atan2(t251, t262);
	t237 = sin(t242);
	t238 = cos(t242);
	t222 = t237 * t251 + t238 * t262;
	t219 = 0.1e1 / t222;
	t229 = 0.1e1 / t293;
	t259 = 0.1e1 / t262;
	t220 = 0.1e1 / t222 ^ 2;
	t260 = 0.1e1 / t262 ^ 2;
	t246 = t251 ^ 2;
	t241 = t246 * t260 + 0.1e1;
	t239 = 0.1e1 / t241;
	t313 = t280 * t282;
	t254 = (t273 * t313 + t274 * t277) * t278;
	t248 = t254 * t275;
	t308 = t273 * t314;
	t320 = t251 * t260;
	t296 = t308 * t320;
	t214 = (-t248 * t259 - t275 * t296) * t239;
	t295 = -t237 * t262 + t238 * t251;
	t307 = t275 * t314;
	t298 = t273 * t307;
	t211 = t295 * t214 - t237 * t248 + t238 * t298;
	t324 = t211 * t219 * t220;
	t317 = t273 * t280;
	t250 = t294 * t277 + t278 * t317;
	t323 = t220 * t250;
	t297 = qJD(5) * t306;
	t299 = -t325 * qJD(5) + t263 * t275;
	t300 = t292 * qJD(5) + t265 * t275;
	t225 = t299 * t273 - t300 * t274 + t284 * t297;
	t231 = t229 * t230;
	t322 = t225 * t231;
	t243 = t263 * t273 - t265 * t274;
	t321 = t235 * t243;
	t316 = t275 * t278;
	t228 = 0.1e1 + t327;
	t301 = t263 * qJD(5) - t275 * t325;
	t224 = t326 * t273 + t301 * t274 - t283 * t297;
	t309 = t235 * t230 * t224;
	t310 = 0.2e1 * (-t232 * t322 + t309) / t228 ^ 2;
	t253 = (-t273 * t315 + t274 * t280) * t278;
	t291 = t254 * t259 + t296;
	t261 = t259 * t260;
	t249 = (t274 * t313 - t318) * t316;
	t247 = t275 * t253;
	t245 = t250 ^ 2;
	t244 = -t273 * t292 + t274 * t325;
	t226 = 0.1e1 / t228;
	t218 = t245 * t220 + 0.1e1;
	t215 = t291 * t239;
	t212 = -t295 * t215 - t237 * t254 + t238 * t308;
	t210 = 0.2e1 * t291 / t241 ^ 2 * (-t246 * t261 * t298 - t248 * t320) + (-t249 * t259 + (0.2e1 * t251 * t261 * t273 ^ 2 * t307 + (-t251 * t274 * t275 + 0.2e1 * t248 * t273) * t260) * t314) * t239;
	t208 = (-t229 * t244 - t230 * t321) * t310 + ((t301 * t273 - t326 * t274) * t229 - 0.2e1 * t321 * t322 + (-t244 * t225 + (t300 * t273 + t299 * t274) * t235 + t243 * t224) * t230) * t226;
	t207 = 0.2e1 * (t212 * t323 - t219 * t253) / t218 ^ 2 * (-t245 * t324 + t247 * t323) + ((-t274 * t315 - t317) * t219 * t316 + (-t253 * t211 - t212 * t247) * t220 + (0.2e1 * t212 * t324 + (-(t275 * t303 + t210 * t251 + t215 * t248 + (t215 * t262 - t254) * t214) * t238 - (t215 * t298 - t210 * t262 - t249 + (t215 * t251 - t308) * t214) * t237) * t220) * t250) / t218;
	t1 = [0, t210, t210, t210, 0; 0, t207, t207, t207, 0; 0, t208, t208, t208, (-t229 * t293 - t327) * t310 + (0.2e1 * t309 + (-t230 * t293 - 0.2e1 * t232 * t231 + t229) * t225) * t226;];
	JaD_rot = t1;
end