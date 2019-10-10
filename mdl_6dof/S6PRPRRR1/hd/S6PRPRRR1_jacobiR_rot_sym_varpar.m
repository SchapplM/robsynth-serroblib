% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:53
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRRR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:47
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:47
	% EndTime: 2019-10-09 21:53:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t22 = cos(pkin(6));
	t23 = sin(qJ(2));
	t26 = t22 * t23;
	t24 = cos(qJ(2));
	t25 = t22 * t24;
	t21 = cos(pkin(11));
	t20 = sin(pkin(6));
	t19 = sin(pkin(11));
	t1 = [0, -t19 * t25 - t21 * t23, 0, 0, 0, 0; 0, -t19 * t23 + t21 * t25, 0, 0, 0, 0; 0, t20 * t24, 0, 0, 0, 0; 0, t19 * t26 - t21 * t24, 0, 0, 0, 0; 0, -t19 * t24 - t21 * t26, 0, 0, 0, 0; 0, -t20 * t23, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (14->7), mult. (40->16), div. (0->0), fcn. (60->8), ass. (0->13)
	t35 = sin(pkin(12));
	t38 = cos(pkin(12));
	t41 = sin(qJ(2));
	t42 = cos(qJ(2));
	t34 = -t42 * t35 - t41 * t38;
	t33 = t41 * t35 - t42 * t38;
	t40 = cos(pkin(6));
	t39 = cos(pkin(11));
	t37 = sin(pkin(6));
	t36 = sin(pkin(11));
	t32 = t34 * t40;
	t31 = t33 * t40;
	t1 = [0, t36 * t31 + t39 * t34, 0, 0, 0, 0; 0, -t39 * t31 + t36 * t34, 0, 0, 0, 0; 0, -t33 * t37, 0, 0, 0, 0; 0, -t36 * t32 + t39 * t33, 0, 0, 0, 0; 0, t39 * t32 + t36 * t33, 0, 0, 0, 0; 0, t34 * t37, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (44->15), mult. (122->36), div. (0->0), fcn. (178->10), ass. (0->23)
	t93 = sin(pkin(6));
	t97 = sin(qJ(4));
	t104 = t93 * t97;
	t99 = cos(qJ(4));
	t103 = t93 * t99;
	t100 = cos(qJ(2));
	t91 = sin(pkin(12));
	t94 = cos(pkin(12));
	t98 = sin(qJ(2));
	t102 = t100 * t94 - t98 * t91;
	t101 = t100 * t91 + t98 * t94;
	t96 = cos(pkin(6));
	t87 = t101 * t96;
	t92 = sin(pkin(11));
	t95 = cos(pkin(11));
	t81 = t102 * t92 + t95 * t87;
	t83 = t102 * t95 - t92 * t87;
	t86 = t102 * t96;
	t85 = t101 * t93;
	t84 = t102 * t93;
	t82 = -t101 * t95 - t92 * t86;
	t80 = -t101 * t92 + t95 * t86;
	t1 = [0, t82 * t99, 0, t103 * t92 - t83 * t97, 0, 0; 0, t80 * t99, 0, -t103 * t95 - t81 * t97, 0, 0; 0, t84 * t99, 0, -t85 * t97 + t96 * t99, 0, 0; 0, -t82 * t97, 0, -t104 * t92 - t83 * t99, 0, 0; 0, -t80 * t97, 0, t104 * t95 - t81 * t99, 0, 0; 0, -t84 * t97, 0, -t85 * t99 - t96 * t97, 0, 0; 0, t83, 0, 0, 0, 0; 0, t81, 0, 0, 0, 0; 0, t85, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (94->16), mult. (178->36), div. (0->0), fcn. (260->10), ass. (0->30)
	t120 = sin(pkin(11));
	t121 = sin(pkin(6));
	t129 = t120 * t121;
	t123 = cos(pkin(11));
	t128 = t121 * t123;
	t124 = cos(pkin(6));
	t119 = sin(pkin(12));
	t122 = cos(pkin(12));
	t125 = sin(qJ(2));
	t126 = cos(qJ(2));
	t127 = t126 * t119 + t125 * t122;
	t112 = t127 * t124;
	t113 = t125 * t119 - t126 * t122;
	t104 = t123 * t112 - t120 * t113;
	t106 = -t120 * t112 - t123 * t113;
	t118 = qJ(4) + qJ(5);
	t117 = cos(t118);
	t116 = sin(t118);
	t111 = t113 * t124;
	t110 = t127 * t121;
	t109 = t113 * t121;
	t108 = -t110 * t117 - t124 * t116;
	t107 = -t110 * t116 + t124 * t117;
	t105 = t120 * t111 - t123 * t127;
	t103 = -t123 * t111 - t120 * t127;
	t102 = -t106 * t117 - t116 * t129;
	t101 = -t106 * t116 + t117 * t129;
	t100 = -t104 * t117 + t116 * t128;
	t99 = -t104 * t116 - t117 * t128;
	t1 = [0, t105 * t117, 0, t101, t101, 0; 0, t103 * t117, 0, t99, t99, 0; 0, -t109 * t117, 0, t107, t107, 0; 0, -t105 * t116, 0, t102, t102, 0; 0, -t103 * t116, 0, t100, t100, 0; 0, t109 * t116, 0, t108, t108, 0; 0, t106, 0, 0, 0, 0; 0, t104, 0, 0, 0, 0; 0, t110, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (204->32), mult. (409->65), div. (0->0), fcn. (583->12), ass. (0->40)
	t180 = sin(pkin(12));
	t183 = cos(pkin(12));
	t187 = sin(qJ(2));
	t189 = cos(qJ(2));
	t173 = t187 * t180 - t189 * t183;
	t179 = qJ(4) + qJ(5);
	t177 = sin(t179);
	t178 = cos(t179);
	t185 = cos(pkin(6));
	t191 = t189 * t180 + t187 * t183;
	t172 = t191 * t185;
	t181 = sin(pkin(11));
	t184 = cos(pkin(11));
	t193 = t184 * t172 - t181 * t173;
	t182 = sin(pkin(6));
	t196 = t182 * t184;
	t156 = -t177 * t193 - t178 * t196;
	t186 = sin(qJ(6));
	t202 = t156 * t186;
	t192 = -t181 * t172 - t184 * t173;
	t197 = t181 * t182;
	t158 = -t177 * t192 + t178 * t197;
	t201 = t158 * t186;
	t171 = t191 * t182;
	t167 = -t171 * t177 + t185 * t178;
	t200 = t167 * t186;
	t199 = t178 * t186;
	t188 = cos(qJ(6));
	t198 = t178 * t188;
	t190 = t173 * t185;
	t170 = t173 * t182;
	t168 = t171 * t178 + t185 * t177;
	t166 = t167 * t188;
	t164 = t181 * t190 - t184 * t191;
	t161 = -t181 * t191 - t184 * t190;
	t159 = t177 * t197 + t178 * t192;
	t157 = -t177 * t196 + t178 * t193;
	t155 = t158 * t188;
	t154 = t156 * t188;
	t1 = [0, t164 * t198 + t186 * t192, 0, t155, t155, -t159 * t186 - t164 * t188; 0, t161 * t198 + t186 * t193, 0, t154, t154, -t157 * t186 - t161 * t188; 0, -t170 * t198 + t171 * t186, 0, t166, t166, -t168 * t186 + t170 * t188; 0, -t164 * t199 + t188 * t192, 0, -t201, -t201, -t159 * t188 + t164 * t186; 0, -t161 * t199 + t188 * t193, 0, -t202, -t202, -t157 * t188 + t161 * t186; 0, t170 * t199 + t171 * t188, 0, -t200, -t200, -t168 * t188 - t170 * t186; 0, t164 * t177, 0, t159, t159, 0; 0, t161 * t177, 0, t157, t157, 0; 0, -t170 * t177, 0, t168, t168, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end