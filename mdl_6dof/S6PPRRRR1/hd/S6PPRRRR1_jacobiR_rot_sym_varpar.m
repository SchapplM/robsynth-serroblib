% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:18
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPRRRR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_jacobiR_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (20->14), mult. (62->33), div. (0->0), fcn. (88->10), ass. (0->20)
	t53 = sin(pkin(12));
	t59 = cos(pkin(6));
	t68 = t53 * t59;
	t54 = sin(pkin(7));
	t55 = sin(pkin(6));
	t67 = t54 * t55;
	t66 = t54 * t59;
	t56 = cos(pkin(13));
	t58 = cos(pkin(7));
	t65 = t56 * t58;
	t57 = cos(pkin(12));
	t64 = t57 * t59;
	t52 = sin(pkin(13));
	t63 = -(-t53 * t52 + t56 * t64) * t58 + t57 * t67;
	t62 = (-t57 * t52 - t56 * t68) * t58 + t53 * t67;
	t61 = cos(qJ(3));
	t60 = sin(qJ(3));
	t51 = -t52 * t68 + t57 * t56;
	t49 = t52 * t64 + t53 * t56;
	t1 = [0, 0, -t51 * t60 + t62 * t61, 0, 0, 0; 0, 0, -t49 * t60 - t63 * t61, 0, 0, 0; 0, 0, t61 * t66 + (-t52 * t60 + t61 * t65) * t55, 0, 0, 0; 0, 0, -t51 * t61 - t62 * t60, 0, 0, 0; 0, 0, -t49 * t61 + t63 * t60, 0, 0, 0; 0, 0, -t60 * t66 + (-t52 * t61 - t60 * t65) * t55, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:21
	% EndTime: 2019-10-09 21:18:21
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (69->26), mult. (203->58), div. (0->0), fcn. (284->12), ass. (0->34)
	t116 = sin(pkin(12));
	t122 = cos(pkin(6));
	t134 = t116 * t122;
	t117 = sin(pkin(7));
	t118 = sin(pkin(6));
	t133 = t117 * t118;
	t132 = t117 * t122;
	t121 = cos(pkin(7));
	t131 = t118 * t121;
	t119 = cos(pkin(13));
	t130 = t119 * t121;
	t120 = cos(pkin(12));
	t129 = t120 * t122;
	t115 = sin(pkin(13));
	t111 = -t116 * t115 + t119 * t129;
	t128 = t111 * t121 - t120 * t133;
	t113 = -t120 * t115 - t119 * t134;
	t127 = t113 * t121 + t116 * t133;
	t126 = cos(qJ(3));
	t125 = cos(qJ(4));
	t124 = sin(qJ(3));
	t123 = sin(qJ(4));
	t114 = -t115 * t134 + t120 * t119;
	t112 = t115 * t129 + t116 * t119;
	t110 = -t119 * t133 + t122 * t121;
	t109 = -t113 * t117 + t116 * t131;
	t108 = -t111 * t117 - t120 * t131;
	t107 = t124 * t132 + (t115 * t126 + t124 * t130) * t118;
	t106 = t126 * t132 + (-t115 * t124 + t126 * t130) * t118;
	t105 = t114 * t126 + t127 * t124;
	t104 = -t114 * t124 + t127 * t126;
	t103 = t112 * t126 + t128 * t124;
	t102 = -t112 * t124 + t128 * t126;
	t1 = [0, 0, t104 * t125, -t105 * t123 + t109 * t125, 0, 0; 0, 0, t102 * t125, -t103 * t123 + t108 * t125, 0, 0; 0, 0, t106 * t125, -t107 * t123 + t110 * t125, 0, 0; 0, 0, -t104 * t123, -t105 * t125 - t109 * t123, 0, 0; 0, 0, -t102 * t123, -t103 * t125 - t108 * t123, 0, 0; 0, 0, -t106 * t123, -t107 * t125 - t110 * t123, 0, 0; 0, 0, t105, 0, 0, 0; 0, 0, t103, 0, 0, 0; 0, 0, t107, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:21
	% EndTime: 2019-10-09 21:18:21
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (135->27), mult. (307->58), div. (0->0), fcn. (430->12), ass. (0->41)
	t144 = sin(pkin(12));
	t150 = cos(pkin(6));
	t160 = t144 * t150;
	t145 = sin(pkin(7));
	t146 = sin(pkin(6));
	t159 = t145 * t146;
	t158 = t145 * t150;
	t149 = cos(pkin(7));
	t157 = t146 * t149;
	t147 = cos(pkin(13));
	t156 = t147 * t149;
	t148 = cos(pkin(12));
	t155 = t148 * t150;
	t143 = sin(pkin(13));
	t136 = -t144 * t143 + t147 * t155;
	t154 = t136 * t149 - t148 * t159;
	t138 = -t148 * t143 - t147 * t160;
	t153 = t138 * t149 + t144 * t159;
	t152 = cos(qJ(3));
	t151 = sin(qJ(3));
	t142 = qJ(4) + qJ(5);
	t141 = cos(t142);
	t140 = sin(t142);
	t139 = -t143 * t160 + t148 * t147;
	t137 = t143 * t155 + t144 * t147;
	t135 = -t147 * t159 + t150 * t149;
	t134 = -t138 * t145 + t144 * t157;
	t133 = -t136 * t145 - t148 * t157;
	t132 = t151 * t158 + (t143 * t152 + t151 * t156) * t146;
	t131 = t152 * t158 + (-t143 * t151 + t152 * t156) * t146;
	t130 = t139 * t152 + t153 * t151;
	t129 = -t139 * t151 + t153 * t152;
	t128 = t137 * t152 + t154 * t151;
	t127 = -t137 * t151 + t154 * t152;
	t126 = -t132 * t141 - t135 * t140;
	t125 = -t132 * t140 + t135 * t141;
	t124 = -t130 * t141 - t134 * t140;
	t123 = -t130 * t140 + t134 * t141;
	t122 = -t128 * t141 - t133 * t140;
	t121 = -t128 * t140 + t133 * t141;
	t1 = [0, 0, t129 * t141, t123, t123, 0; 0, 0, t127 * t141, t121, t121, 0; 0, 0, t131 * t141, t125, t125, 0; 0, 0, -t129 * t140, t124, t124, 0; 0, 0, -t127 * t140, t122, t122, 0; 0, 0, -t131 * t140, t126, t126, 0; 0, 0, t130, 0, 0, 0; 0, 0, t128, 0, 0, 0; 0, 0, t132, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:21
	% EndTime: 2019-10-09 21:18:21
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (289->47), mult. (678->91), div. (0->0), fcn. (937->14), ass. (0->53)
	t239 = cos(qJ(3));
	t211 = sin(pkin(13));
	t212 = sin(pkin(12));
	t215 = cos(pkin(13));
	t216 = cos(pkin(12));
	t218 = cos(pkin(6));
	t228 = t216 * t218;
	t203 = t211 * t228 + t212 * t215;
	t220 = sin(qJ(3));
	t217 = cos(pkin(7));
	t225 = -t212 * t211 + t215 * t228;
	t223 = t225 * t217;
	t213 = sin(pkin(7));
	t214 = sin(pkin(6));
	t231 = t214 * t213;
	t193 = t203 * t239 + (-t216 * t231 + t223) * t220;
	t230 = t214 * t217;
	t198 = -t225 * t213 - t216 * t230;
	t210 = qJ(4) + qJ(5);
	t208 = sin(t210);
	t209 = cos(t210);
	t185 = -t193 * t208 + t198 * t209;
	t219 = sin(qJ(6));
	t238 = t185 * t219;
	t233 = t212 * t218;
	t204 = -t211 * t233 + t216 * t215;
	t224 = t216 * t211 + t215 * t233;
	t222 = t224 * t217;
	t195 = t204 * t239 + (t212 * t231 - t222) * t220;
	t199 = t212 * t230 + t224 * t213;
	t187 = -t195 * t208 + t199 * t209;
	t237 = t187 * t219;
	t229 = t215 * t217;
	t232 = t213 * t218;
	t197 = t220 * t232 + (t239 * t211 + t220 * t229) * t214;
	t202 = -t215 * t231 + t218 * t217;
	t190 = -t197 * t208 + t202 * t209;
	t236 = t190 * t219;
	t235 = t209 * t219;
	t221 = cos(qJ(6));
	t234 = t209 * t221;
	t227 = t214 * t239;
	t226 = t213 * t227;
	t196 = t214 * t211 * t220 - t227 * t229 - t239 * t232;
	t194 = t204 * t220 - t212 * t226 + t239 * t222;
	t192 = t203 * t220 + t216 * t226 - t239 * t223;
	t191 = t197 * t209 + t202 * t208;
	t189 = t190 * t221;
	t188 = t195 * t209 + t199 * t208;
	t186 = t193 * t209 + t198 * t208;
	t184 = t187 * t221;
	t183 = t185 * t221;
	t1 = [0, 0, -t194 * t234 + t195 * t219, t184, t184, -t188 * t219 + t194 * t221; 0, 0, -t192 * t234 + t193 * t219, t183, t183, -t186 * t219 + t192 * t221; 0, 0, -t196 * t234 + t197 * t219, t189, t189, -t191 * t219 + t196 * t221; 0, 0, t194 * t235 + t195 * t221, -t237, -t237, -t188 * t221 - t194 * t219; 0, 0, t192 * t235 + t193 * t221, -t238, -t238, -t186 * t221 - t192 * t219; 0, 0, t196 * t235 + t197 * t221, -t236, -t236, -t191 * t221 - t196 * t219; 0, 0, -t194 * t208, t188, t188, 0; 0, 0, -t192 * t208, t186, t186, 0; 0, 0, -t196 * t208, t191, t191, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end