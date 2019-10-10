% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6RPRPRR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:53
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:11
	% DurationCPUTime: 0.81s
	% Computational Cost: add. (1893->71), mult. (1835->159), div. (470->13), fcn. (2177->7), ass. (0->71)
	t89 = sin(qJ(1));
	t119 = qJD(1) * t89;
	t138 = 0.2e1 * t89;
	t116 = qJD(3) * t89;
	t90 = cos(qJ(1));
	t118 = qJD(1) * t90;
	t82 = pkin(10) + qJ(3);
	t80 = sin(t82);
	t109 = t80 * t118;
	t76 = t80 ^ 2;
	t81 = cos(t82);
	t78 = 0.1e1 / t81 ^ 2;
	t124 = t76 * t78;
	t84 = t89 ^ 2;
	t71 = t84 * t124 + 0.1e1;
	t69 = 0.1e1 / t71;
	t77 = 0.1e1 / t81;
	t55 = (-(-t81 * t116 - t109) * t77 + t116 * t124) * t69;
	t136 = t55 - t116;
	t85 = t90 ^ 2;
	t86 = 0.1e1 / t90;
	t135 = (t84 / t85 + 0.1e1) * t86 * t119;
	t120 = t89 * t80;
	t68 = atan2(-t120, -t81);
	t66 = sin(t68);
	t111 = t66 * t120;
	t67 = cos(t68);
	t62 = -t67 * t81 - t111;
	t59 = 0.1e1 / t62;
	t60 = 0.1e1 / t62 ^ 2;
	t134 = -0.2e1 * t80;
	t133 = t69 - 0.1e1;
	t126 = t67 * t80;
	t51 = (-t55 * t89 + qJD(3)) * t126 + (t136 * t81 - t109) * t66;
	t132 = t51 * t59 * t60;
	t131 = t55 * t80;
	t130 = t60 * t80;
	t129 = t60 * t90;
	t102 = t76 * t89 * t118;
	t122 = t77 * t80;
	t75 = t80 * t76;
	t79 = t77 * t78;
	t98 = qJD(3) * (t75 * t79 + t122);
	t128 = (t78 * t102 + t84 * t98) / t71 ^ 2;
	t127 = t66 * t89;
	t125 = t76 * t77;
	t123 = t76 * t85;
	t121 = t84 / t90 ^ 2;
	t117 = qJD(3) * t81;
	t58 = t60 * t123 + 0.1e1;
	t115 = 0.2e1 * (-t123 * t132 + (t80 * t85 * t117 - t102) * t60) / t58 ^ 2;
	t114 = 0.2e1 * t132;
	t74 = t78 * t121 + 0.1e1;
	t113 = 0.2e1 * (t79 * qJD(3) * t80 * t121 + t78 * t135) / t74 ^ 2;
	t112 = t80 * t129;
	t110 = t69 * t125;
	t108 = 0.1e1 + t124;
	t107 = 0.1e1 + t121;
	t106 = t80 * t115;
	t105 = t128 * t134;
	t104 = t128 * t138;
	t103 = t89 * t110;
	t101 = t108 * t90;
	t99 = t107 * t80 * t78;
	t72 = 0.1e1 / t74;
	t64 = t108 * t89 * t69;
	t56 = 0.1e1 / t58;
	t54 = (t133 * t80 * t66 - t67 * t103) * t90;
	t53 = -t81 * t127 + t126 + (-t67 * t120 + t66 * t81) * t64;
	t52 = -t108 * t104 + (qJD(1) * t101 + t98 * t138) * t69;
	t1 = [t90 * t77 * t105 + (qJD(3) * t101 - t119 * t122) * t69, 0, t52, 0, 0, 0; (t59 * t106 + (-t59 * t117 + (qJD(1) * t54 + t51) * t130) * t56) * t89 + (t60 * t106 * t54 + (-((t55 * t103 + t133 * t117 + t105) * t66 + (t104 * t125 - t131 + (t131 + (-t75 * t78 + t134) * t116) * t69) * t67) * t112 + (t80 * t114 - t60 * t117) * t54 + (-t59 + ((-t84 + t85) * t67 * t110 + t133 * t111) * t60) * t80 * qJD(1)) * t56) * t90, 0, (t53 * t130 - t59 * t81) * t90 * t115 + ((-t59 * t119 + (-qJD(3) * t53 - t51) * t129) * t81 + (-t90 * qJD(3) * t59 - (-t52 * t67 * t89 - t136 * t66 + (-qJD(3) * t66 - t118 * t67 + t127 * t55) * t64) * t112 + (t90 * t114 + t60 * t119) * t53 - ((t52 - t118) * t66 + ((-t64 * t89 + 0.1e1) * qJD(3) + (t64 - t89) * t55) * t67) * t81 * t129) * t80) * t56, 0, 0, 0; t107 * t77 * t113 + (-qJD(3) * t99 - 0.2e1 * t77 * t135) * t72, 0, t86 * t78 * t113 * t120 + ((-0.2e1 * t76 * t79 - t77) * t86 * t116 - qJD(1) * t99) * t72, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:12
	% DurationCPUTime: 1.78s
	% Computational Cost: add. (7491->94), mult. (9702->187), div. (711->12), fcn. (12089->11), ass. (0->95)
	t271 = pkin(10) + qJ(3);
	t260 = sin(t271);
	t261 = cos(t271);
	t290 = sin(qJ(5));
	t292 = cos(qJ(5));
	t208 = t260 * t292 - t261 * t290;
	t291 = sin(qJ(1));
	t199 = t208 * t291;
	t207 = t260 * t290 + t261 * t292;
	t186 = atan2(t199, t207);
	t181 = sin(t186);
	t182 = cos(t186);
	t201 = t207 * t291;
	t255 = -t181 * t207 + t182 * t199;
	t197 = t199 ^ 2;
	t205 = 0.1e1 / t207 ^ 2;
	t185 = t197 * t205 + 0.1e1;
	t183 = 0.1e1 / t185;
	t204 = 0.1e1 / t207;
	t273 = t199 * t208;
	t247 = -t201 * t204 - t205 * t273;
	t300 = t247 * t183;
	t158 = -t181 * t201 + t182 * t208 + t255 * t300;
	t293 = cos(qJ(1));
	t203 = t208 * t293;
	t319 = t158 * t203;
	t171 = t181 * t199 + t182 * t207;
	t168 = 0.1e1 / t171;
	t169 = 0.1e1 / t171 ^ 2;
	t202 = t207 * t293;
	t198 = t203 ^ 2;
	t167 = t198 * t169 + 0.1e1;
	t308 = qJD(3) - qJD(5);
	t175 = -t199 * qJD(1) + t308 * t202;
	t281 = t175 * t169;
	t176 = -t203 * qJD(1) - t308 * t201;
	t188 = t308 * t208;
	t274 = t199 * t205;
	t249 = -t176 * t204 + t188 * t274;
	t161 = t249 * t183;
	t156 = t255 * t161 - t181 * t176 - t182 * t188;
	t289 = t156 * t168 * t169;
	t270 = 0.2e1 * (-t198 * t289 + t203 * t281) / t167 ^ 2;
	t318 = (t168 * t202 + t169 * t319) * t270;
	t174 = qJD(1) * t201 + t308 * t203;
	t226 = sin(qJ(6));
	t227 = cos(qJ(6));
	t196 = t202 * t227 - t291 * t226;
	t263 = qJD(1) * t293;
	t172 = t196 * qJD(6) - t174 * t226 + t227 * t263;
	t245 = -t202 * t226 - t291 * t227;
	t309 = t245 * qJD(6);
	t173 = -t174 * t227 - t226 * t263 + t309;
	t189 = t245 ^ 2;
	t191 = 0.1e1 / t196 ^ 2;
	t180 = t189 * t191 + 0.1e1;
	t178 = 0.1e1 / t180;
	t190 = 0.1e1 / t196;
	t276 = t191 * t245;
	t248 = -t226 * t190 - t227 * t276;
	t283 = t173 * t190 * t191;
	t295 = -0.2e1 * t245;
	t259 = t283 * t295;
	t317 = (t248 * t175 - (((-t173 - t309) * t226 - t172 * t227) * t191 + (qJD(6) * t190 + t259) * t227) * t203) * t178;
	t316 = -t202 * t156 + t158 * t175;
	t269 = 0.2e1 * t289;
	t315 = -t174 * t168 - t269 * t319;
	t187 = t308 * t207;
	t313 = (t207 * t300 + t201) * t161 + t300 * t176 - t187;
	t177 = qJD(1) * t202 - t308 * t199;
	t312 = -(-t199 * t300 - t208) * t161 - t300 * t188 + t177;
	t302 = t188 * t205;
	t277 = t204 * t302;
	t310 = ((t176 * t208 - t187 * t199 - t188 * t201) * t205 - t177 * t204 - 0.2e1 * t273 * t277) * t183;
	t275 = t199 * t204;
	t298 = (-t182 * t275 + t181) * t183 - t181;
	t294 = -0.2e1 * t203;
	t288 = (-t172 * t276 - t189 * t283) / t180 ^ 2;
	t287 = (-t176 * t274 + t197 * t277) / t185 ^ 2;
	t280 = t178 * t191;
	t279 = t181 * t203;
	t278 = t182 * t203;
	t268 = -0.2e1 * t288;
	t267 = -0.2e1 * t287;
	t266 = t191 * t288;
	t265 = t204 * t287;
	t264 = t172 * t280;
	t262 = qJD(1) * t291;
	t194 = -t201 * t227 - t293 * t226;
	t246 = t201 * t226 - t293 * t227;
	t165 = 0.1e1 / t167;
	t160 = t298 * t203;
	t154 = t247 * t267 + t310;
	t153 = 0.2e1 * t247 * t287 - t310;
	t1 = [t265 * t294 + (t175 * t204 + t203 * t302) * t183, 0, t153, 0, t154, 0; -t199 * t168 * t270 + (-t176 * t168 + (-t156 * t199 - t160 * t175) * t169) * t165 - ((-t160 * t269 + t298 * t281) * t165 + (-t160 * t270 + ((t161 * t183 * t275 + t267) * t279 + (0.2e1 * t199 * t265 - t161 + (t161 - t249) * t183) * t278) * t165) * t169) * t203, 0, t318 + (((t153 * t199 + t313) * t278 + (-t153 * t207 + t312) * t279 - t316) * t169 - t315) * t165, 0, -t318 + (((t154 * t199 - t313) * t278 + (-t154 * t207 - t312) * t279 + t316) * t169 + t315) * t165, 0; (t266 * t295 - t264) * t194 - (-t173 * t280 + t190 * t268) * t246 + ((t194 * qJD(6) - t177 * t226 - t227 * t262) * t190 + (t246 * qJD(6) - t177 * t227 + t226 * t262) * t276 + t194 * t259) * t178, 0, t248 * t288 * t294 + t317, 0, -t248 * t203 * t268 - t317, t268 + (t264 - (-t178 * t283 - t266) * t245) * t295;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end