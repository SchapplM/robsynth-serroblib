% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPPPR2
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
%   Wie in S5RPPPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPPPR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
JaD_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (265->32), mult. (614->97), div. (108->12), fcn. (792->9), ass. (0->49)
	t86 = sin(pkin(7));
	t79 = t86 ^ 2;
	t88 = cos(pkin(7));
	t81 = 0.1e1 / t88 ^ 2;
	t89 = sin(qJ(1));
	t83 = t89 ^ 2;
	t77 = t79 * t81 * t83 + 0.1e1;
	t90 = cos(qJ(1));
	t84 = t90 ^ 2;
	t106 = 0.1e1 / t77 ^ 2 * t84;
	t112 = t106 * t81;
	t75 = 0.1e1 / t77;
	t111 = (t75 - 0.1e1) * t86;
	t104 = t86 * t89;
	t74 = atan2(-t104, -t88);
	t72 = sin(t74);
	t73 = cos(t74);
	t58 = -t72 * t104 - t73 * t88;
	t55 = 0.1e1 / t58;
	t85 = sin(pkin(8));
	t102 = t89 * t85;
	t103 = t88 * t90;
	t87 = cos(pkin(8));
	t71 = t87 * t103 + t102;
	t65 = 0.1e1 / t71;
	t80 = 0.1e1 / t88;
	t56 = 0.1e1 / t58 ^ 2;
	t66 = 0.1e1 / t71 ^ 2;
	t109 = t56 * t90;
	t101 = t89 * t87;
	t70 = t85 * t103 - t101;
	t108 = t66 * t70;
	t69 = -t88 * t101 + t85 * t90;
	t107 = t69 * t70;
	t105 = t79 * t80;
	t100 = qJD(1) * t89;
	t99 = t80 * t112;
	t68 = -t88 * t102 - t87 * t90;
	t51 = (-t73 * t75 * t89 * t105 + t72 * t111) * t90;
	t78 = t86 * t79;
	t67 = t65 * t66;
	t64 = t70 ^ 2;
	t63 = t69 * qJD(1);
	t62 = t68 * qJD(1);
	t61 = t64 * t66 + 0.1e1;
	t57 = t55 * t56;
	t54 = t56 * t79 * t84 + 0.1e1;
	t50 = qJD(1) * t51;
	t1 = [(-t75 * t80 * t86 - 0.2e1 * t78 * t99) * t100, 0, 0, 0, 0; (0.2e1 * (t51 * t109 + t55 * t89) / t54 ^ 2 * (-t50 * t57 * t84 - t100 * t109) * t79 + ((0.2e1 * t51 * t57 * t90 + t56 * t89) * t50 + (-t90 * t55 + ((t51 + (t78 * t112 + t111) * t90 * t72) * t89 - (0.2e1 * t83 * t79 ^ 2 * t99 + (t106 + (t83 - 0.2e1 * t84) * t75) * t105) * t90 * t73) * t56) * qJD(1)) / t54) * t86, 0, 0, 0, 0; 0.2e1 * (t66 * t107 - t65 * t68) / t61 ^ 2 * (-t63 * t64 * t67 + t62 * t108) + (-t69 * t62 * t66 + (0.2e1 * t67 * t107 - t68 * t66) * t63 + (t71 * t108 - t70 * t65) * qJD(1)) / t61, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:32
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (367->38), mult. (1253->100), div. (154->13), fcn. (1694->11), ass. (0->60)
	t122 = sin(pkin(8));
	t123 = sin(pkin(7));
	t149 = t123 * t122;
	t151 = 0.1e1 / t149;
	t150 = 0.1e1 / t122 ^ 2 / t123 ^ 2;
	t125 = cos(pkin(7));
	t153 = cos(pkin(8));
	t158 = sin(qJ(1));
	t137 = t158 * t153;
	t126 = cos(qJ(1));
	t147 = t126 * t122;
	t111 = t125 * t147 - t137;
	t139 = t126 * t153;
	t141 = t158 * t122;
	t108 = t125 * t141 + t139;
	t101 = atan2(-t108, t149);
	t98 = cos(t101);
	t146 = t98 * t151;
	t97 = sin(t101);
	t102 = t108 ^ 2 * t150 + 0.1e1;
	t99 = 0.1e1 / t102;
	t135 = (t108 * t146 + t97) * t99 - t97;
	t83 = -t97 * t108 + t98 * t149;
	t80 = 0.1e1 / t83;
	t112 = t125 * t139 + t141;
	t121 = sin(pkin(9));
	t124 = cos(pkin(9));
	t148 = t123 * t126;
	t96 = t112 * t124 + t121 * t148;
	t90 = 0.1e1 / t96;
	t81 = 0.1e1 / t83 ^ 2;
	t91 = 0.1e1 / t96 ^ 2;
	t105 = t111 * qJD(1);
	t74 = t135 * t105;
	t157 = t74 * t80 * t81;
	t110 = -t125 * t137 + t147;
	t104 = t110 * qJD(1);
	t142 = t123 * t158;
	t138 = qJD(1) * t142;
	t88 = t104 * t124 - t121 * t138;
	t156 = t88 * t90 * t91;
	t94 = t110 * t124 - t121 * t142;
	t95 = t112 * t121 - t124 * t148;
	t155 = t94 * t95;
	t154 = t111 * t81;
	t152 = t105 * t111;
	t145 = t99 * t151;
	t144 = t108 * t150 * t151;
	t140 = qJD(1) * t148;
	t107 = t111 ^ 2;
	t106 = t112 * qJD(1);
	t103 = t108 * qJD(1);
	t100 = 0.1e1 / t102 ^ 2;
	t93 = t110 * t121 + t124 * t142;
	t89 = t95 ^ 2;
	t87 = t104 * t121 + t124 * t138;
	t86 = t89 * t91 + 0.1e1;
	t79 = t107 * t81 + 0.1e1;
	t75 = t135 * t111;
	t1 = [0.2e1 * t100 * t144 * t152 + t103 * t145, 0, 0, 0, 0; 0.2e1 * (t108 * t80 + t75 * t154) / t79 ^ 2 * (-t103 * t154 - t107 * t157) + (0.2e1 * t75 * t157 * t111 - t105 * t80 + (t75 * t103 + t108 * t74 + (t135 * t103 - (0.2e1 * t98 * t145 + (-t146 + (-0.2e1 * t98 * t144 - t97 * t150) * t108) * t100) * t152) * t111) * t81) / t79, 0, 0, 0, 0; 0.2e1 * (t91 * t155 - t90 * t93) / t86 ^ 2 * (t95 * t91 * t87 - t89 * t156) + ((-t106 * t121 + t124 * t140) * t90 + 0.2e1 * t155 * t156 + (-t93 * t88 - (-t106 * t124 - t121 * t140) * t95 - t94 * t87) * t91) / t86, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:32
	% EndTime: 2022-01-23 09:00:32
	% DurationCPUTime: 0.49s
	% Computational Cost: add. (1052->54), mult. (3140->129), div. (139->12), fcn. (4132->13), ass. (0->72)
	t201 = cos(pkin(8));
	t199 = sin(pkin(9));
	t245 = sin(pkin(7));
	t221 = t245 * t199;
	t202 = cos(pkin(7));
	t246 = cos(pkin(9));
	t222 = t202 * t246;
	t191 = t201 * t222 + t221;
	t203 = sin(qJ(5));
	t205 = cos(qJ(5));
	t200 = sin(pkin(8));
	t230 = t200 * t202;
	t174 = t191 * t205 + t203 * t230;
	t223 = t200 * t246;
	t193 = -t203 * t201 + t205 * t223;
	t204 = sin(qJ(1));
	t206 = cos(qJ(1));
	t157 = -t174 * t204 + t193 * t206;
	t192 = t205 * t201 + t203 * t223;
	t251 = -t191 * t203 + t205 * t230;
	t159 = t192 * t206 + t204 * t251;
	t171 = t174 * qJD(5);
	t185 = t193 * qJD(5);
	t148 = t159 * qJD(1) + t171 * t206 + t204 * t185;
	t255 = t174 * t206 + t204 * t193;
	t153 = 0.1e1 / t255 ^ 2;
	t248 = -t204 * t192 + t206 * t251;
	t239 = t153 * t248;
	t256 = t148 * t239;
	t218 = t245 * t246;
	t229 = t201 * t202;
	t180 = -t206 * t218 + (t204 * t200 + t206 * t229) * t199;
	t170 = t251 * qJD(5);
	t184 = t192 * qJD(5);
	t147 = t157 * qJD(1) + t170 * t206 - t204 * t184;
	t254 = t147 * t153;
	t177 = (-t200 * t206 + t204 * t229) * t199 - t204 * t218;
	t190 = t201 * t221 + t222;
	t167 = atan2(-t177, t190);
	t162 = sin(t167);
	t163 = cos(t167);
	t146 = -t162 * t177 + t163 * t190;
	t143 = 0.1e1 / t146;
	t152 = 0.1e1 / t255;
	t247 = t177 ^ 2;
	t187 = 0.1e1 / t190;
	t144 = 0.1e1 / t146 ^ 2;
	t188 = 0.1e1 / t190 ^ 2;
	t155 = t248 ^ 2;
	t151 = t153 * t155 + 0.1e1;
	t241 = t152 * t254;
	t244 = (-t155 * t241 - t256) / t151 ^ 2;
	t169 = t180 * qJD(1);
	t166 = t247 * t188 + 0.1e1;
	t164 = 0.1e1 / t166;
	t213 = -t162 + (t163 * t177 * t187 + t162) * t164;
	t138 = t213 * t169;
	t243 = t138 * t143 * t144;
	t242 = t144 * t180;
	t238 = t164 * t187;
	t165 = 0.1e1 / t166 ^ 2;
	t237 = t165 * t177;
	t236 = t169 * t180;
	t224 = 0.2e1 * t244;
	t220 = -0.2e1 * t248 * t241;
	t189 = t187 * t188;
	t172 = t180 ^ 2;
	t168 = t177 * qJD(1);
	t149 = 0.1e1 / t151;
	t142 = t144 * t172 + 0.1e1;
	t139 = t213 * t180;
	t1 = [0.2e1 * t189 * t236 * t237 + t168 * t238, 0, 0, 0, 0; 0.2e1 * (t139 * t242 + t143 * t177) / t142 ^ 2 * (-t168 * t242 - t172 * t243) + (-t169 * t143 + (t138 * t177 + t139 * t168) * t144 + (0.2e1 * t139 * t243 + (t213 * t168 - (-t162 * t188 * t237 + (0.2e1 * t238 + (-0.2e1 * t247 * t189 - t187) * t165) * t163) * t236) * t144) * t180) / t142, 0, 0, 0, 0; (-t152 * t159 - t157 * t239) * t224 + ((qJD(1) * t248 - t171 * t204 + t185 * t206) * t152 + t157 * t220 + (-t159 * t147 + (-qJD(1) * t255 - t170 * t204 - t184 * t206) * t248 - t157 * t148) * t153) * t149, 0, 0, 0, -0.2e1 * t255 * t152 * t244 - t248 * t224 * t239 + ((-t148 * t153 + t220) * t248 + t147 * t152 - t255 * t254 - t256) * t149;];
	JaD_rot = t1;
end