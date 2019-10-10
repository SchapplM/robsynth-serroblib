% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:33
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRPR3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t22 = cos(pkin(6));
	t23 = sin(qJ(2));
	t26 = t22 * t23;
	t24 = cos(qJ(2));
	t25 = t22 * t24;
	t21 = cos(pkin(10));
	t20 = sin(pkin(6));
	t19 = sin(pkin(10));
	t1 = [0, -t19 * t25 - t21 * t23, 0, 0, 0, 0; 0, -t19 * t23 + t21 * t25, 0, 0, 0, 0; 0, t20 * t24, 0, 0, 0, 0; 0, t19 * t26 - t21 * t24, 0, 0, 0, 0; 0, -t19 * t24 - t21 * t26, 0, 0, 0, 0; 0, -t20 * t23, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (14->7), mult. (40->16), div. (0->0), fcn. (60->8), ass. (0->13)
	t35 = sin(pkin(11));
	t38 = cos(pkin(11));
	t41 = sin(qJ(2));
	t42 = cos(qJ(2));
	t34 = -t42 * t35 - t41 * t38;
	t33 = t41 * t35 - t42 * t38;
	t40 = cos(pkin(6));
	t39 = cos(pkin(10));
	t37 = sin(pkin(6));
	t36 = sin(pkin(10));
	t32 = t34 * t40;
	t31 = t33 * t40;
	t1 = [0, t36 * t31 + t39 * t34, 0, 0, 0, 0; 0, -t39 * t31 + t36 * t34, 0, 0, 0, 0; 0, -t33 * t37, 0, 0, 0, 0; 0, -t36 * t32 + t39 * t33, 0, 0, 0, 0; 0, t39 * t32 + t36 * t33, 0, 0, 0, 0; 0, t34 * t37, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (44->15), mult. (122->36), div. (0->0), fcn. (178->10), ass. (0->23)
	t93 = sin(pkin(6));
	t97 = sin(qJ(4));
	t104 = t93 * t97;
	t99 = cos(qJ(4));
	t103 = t93 * t99;
	t100 = cos(qJ(2));
	t91 = sin(pkin(11));
	t94 = cos(pkin(11));
	t98 = sin(qJ(2));
	t102 = t100 * t94 - t98 * t91;
	t101 = t100 * t91 + t98 * t94;
	t96 = cos(pkin(6));
	t87 = t101 * t96;
	t92 = sin(pkin(10));
	t95 = cos(pkin(10));
	t81 = t102 * t92 + t95 * t87;
	t83 = t102 * t95 - t92 * t87;
	t86 = t102 * t96;
	t85 = t101 * t93;
	t84 = t102 * t93;
	t82 = -t101 * t95 - t92 * t86;
	t80 = -t101 * t92 + t95 * t86;
	t1 = [0, t82 * t99, 0, t92 * t103 - t83 * t97, 0, 0; 0, t80 * t99, 0, -t95 * t103 - t81 * t97, 0, 0; 0, t84 * t99, 0, -t85 * t97 + t96 * t99, 0, 0; 0, -t82 * t97, 0, -t92 * t104 - t83 * t99, 0, 0; 0, -t80 * t97, 0, t95 * t104 - t81 * t99, 0, 0; 0, -t84 * t97, 0, -t85 * t99 - t96 * t97, 0, 0; 0, t83, 0, 0, 0, 0; 0, t81, 0, 0, 0, 0; 0, t85, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (44->15), mult. (122->36), div. (0->0), fcn. (178->10), ass. (0->23)
	t114 = sin(pkin(6));
	t118 = sin(qJ(4));
	t124 = t114 * t118;
	t120 = cos(qJ(4));
	t123 = t114 * t120;
	t117 = cos(pkin(6));
	t112 = sin(pkin(11));
	t115 = cos(pkin(11));
	t119 = sin(qJ(2));
	t121 = cos(qJ(2));
	t122 = t121 * t112 + t119 * t115;
	t108 = t122 * t117;
	t109 = t119 * t112 - t121 * t115;
	t113 = sin(pkin(10));
	t116 = cos(pkin(10));
	t102 = t116 * t108 - t113 * t109;
	t104 = -t113 * t108 - t116 * t109;
	t107 = t109 * t117;
	t106 = t122 * t114;
	t105 = t109 * t114;
	t103 = t113 * t107 - t116 * t122;
	t101 = -t116 * t107 - t113 * t122;
	t1 = [0, t104, 0, 0, 0, 0; 0, t102, 0, 0, 0, 0; 0, t106, 0, 0, 0, 0; 0, -t103 * t120, 0, t104 * t118 - t113 * t123, 0, 0; 0, -t101 * t120, 0, t102 * t118 + t116 * t123, 0, 0; 0, t105 * t120, 0, t106 * t118 - t117 * t120, 0, 0; 0, t103 * t118, 0, t104 * t120 + t113 * t124, 0, 0; 0, t101 * t118, 0, t102 * t120 - t116 * t124, 0, 0; 0, -t105 * t118, 0, t106 * t120 + t117 * t118, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:41
	% EndTime: 2019-10-09 21:33:41
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (111->28), mult. (319->65), div. (0->0), fcn. (454->12), ass. (0->33)
	t147 = sin(pkin(11));
	t150 = cos(pkin(11));
	t155 = sin(qJ(2));
	t158 = cos(qJ(2));
	t143 = t155 * t147 - t158 * t150;
	t149 = sin(pkin(6));
	t154 = sin(qJ(4));
	t168 = t149 * t154;
	t157 = cos(qJ(4));
	t167 = t149 * t157;
	t153 = sin(qJ(6));
	t166 = t153 * t154;
	t156 = cos(qJ(6));
	t165 = t154 * t156;
	t152 = cos(pkin(6));
	t160 = t158 * t147 + t155 * t150;
	t142 = t160 * t152;
	t148 = sin(pkin(10));
	t151 = cos(pkin(10));
	t162 = t151 * t142 - t148 * t143;
	t161 = -t148 * t142 - t151 * t143;
	t159 = t143 * t152;
	t141 = t160 * t149;
	t140 = t143 * t149;
	t138 = t141 * t157 + t152 * t154;
	t137 = t141 * t154 - t152 * t157;
	t135 = t148 * t159 - t151 * t160;
	t132 = -t148 * t160 - t151 * t159;
	t130 = t148 * t168 + t157 * t161;
	t129 = -t148 * t167 + t154 * t161;
	t128 = -t151 * t168 + t157 * t162;
	t127 = t151 * t167 + t154 * t162;
	t1 = [0, t135 * t166 + t156 * t161, 0, t130 * t153, 0, t129 * t156 + t135 * t153; 0, t132 * t166 + t156 * t162, 0, t128 * t153, 0, t127 * t156 + t132 * t153; 0, -t140 * t166 + t141 * t156, 0, t138 * t153, 0, t137 * t156 - t140 * t153; 0, t135 * t165 - t153 * t161, 0, t130 * t156, 0, -t129 * t153 + t135 * t156; 0, t132 * t165 - t153 * t162, 0, t128 * t156, 0, -t127 * t153 + t132 * t156; 0, -t140 * t165 - t141 * t153, 0, t138 * t156, 0, -t137 * t153 - t140 * t156; 0, t135 * t157, 0, -t129, 0, 0; 0, t132 * t157, 0, -t127, 0, 0; 0, -t140 * t157, 0, -t137, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end