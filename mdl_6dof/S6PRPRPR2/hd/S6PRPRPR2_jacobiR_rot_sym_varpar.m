% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:31
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRPR2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
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
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.06s
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
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.12s
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
	t1 = [0, t82 * t99, 0, t103 * t92 - t83 * t97, 0, 0; 0, t80 * t99, 0, -t103 * t95 - t81 * t97, 0, 0; 0, t84 * t99, 0, -t85 * t97 + t96 * t99, 0, 0; 0, -t82 * t97, 0, -t104 * t92 - t83 * t99, 0, 0; 0, -t80 * t97, 0, t104 * t95 - t81 * t99, 0, 0; 0, -t84 * t97, 0, -t85 * t99 - t96 * t97, 0, 0; 0, t83, 0, 0, 0, 0; 0, t81, 0, 0, 0, 0; 0, t85, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (74->22), mult. (211->53), div. (0->0), fcn. (300->12), ass. (0->30)
	t130 = sin(pkin(12));
	t140 = cos(qJ(4));
	t148 = t130 * t140;
	t133 = sin(pkin(6));
	t138 = sin(qJ(4));
	t147 = t133 * t138;
	t146 = t133 * t140;
	t134 = cos(pkin(12));
	t145 = t134 * t140;
	t137 = cos(pkin(6));
	t131 = sin(pkin(11));
	t135 = cos(pkin(11));
	t139 = sin(qJ(2));
	t141 = cos(qJ(2));
	t142 = t141 * t131 + t139 * t135;
	t126 = t142 * t137;
	t127 = t139 * t131 - t141 * t135;
	t132 = sin(pkin(10));
	t136 = cos(pkin(10));
	t144 = t136 * t126 - t132 * t127;
	t143 = -t132 * t126 - t136 * t127;
	t125 = t127 * t137;
	t124 = t142 * t133;
	t123 = t127 * t133;
	t122 = -t124 * t138 + t137 * t140;
	t120 = t132 * t125 - t136 * t142;
	t117 = -t136 * t125 - t132 * t142;
	t115 = t132 * t146 - t138 * t143;
	t114 = -t136 * t146 - t138 * t144;
	t1 = [0, t120 * t145 + t130 * t143, 0, t115 * t134, 0, 0; 0, t117 * t145 + t130 * t144, 0, t114 * t134, 0, 0; 0, -t123 * t145 + t124 * t130, 0, t122 * t134, 0, 0; 0, -t120 * t148 + t134 * t143, 0, -t115 * t130, 0, 0; 0, -t117 * t148 + t134 * t144, 0, -t114 * t130, 0, 0; 0, t123 * t148 + t124 * t134, 0, -t122 * t130, 0, 0; 0, t120 * t138, 0, t132 * t147 + t140 * t143, 0, 0; 0, t117 * t138, 0, -t136 * t147 + t140 * t144, 0, 0; 0, -t123 * t138, 0, t124 * t140 + t137 * t138, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (144->29), mult. (319->65), div. (0->0), fcn. (454->12), ass. (0->34)
	t155 = sin(pkin(11));
	t158 = cos(pkin(11));
	t162 = sin(qJ(2));
	t164 = cos(qJ(2));
	t148 = t162 * t155 - t164 * t158;
	t154 = pkin(12) + qJ(6);
	t152 = sin(t154);
	t163 = cos(qJ(4));
	t174 = t152 * t163;
	t153 = cos(t154);
	t173 = t153 * t163;
	t157 = sin(pkin(6));
	t161 = sin(qJ(4));
	t172 = t157 * t161;
	t171 = t157 * t163;
	t160 = cos(pkin(6));
	t166 = t164 * t155 + t162 * t158;
	t147 = t166 * t160;
	t156 = sin(pkin(10));
	t159 = cos(pkin(10));
	t168 = t159 * t147 - t156 * t148;
	t167 = -t156 * t147 - t159 * t148;
	t165 = t148 * t160;
	t146 = t166 * t157;
	t145 = t148 * t157;
	t143 = t146 * t163 + t160 * t161;
	t142 = -t146 * t161 + t160 * t163;
	t140 = t156 * t165 - t159 * t166;
	t137 = -t156 * t166 - t159 * t165;
	t135 = t156 * t172 + t163 * t167;
	t134 = t156 * t171 - t161 * t167;
	t133 = -t159 * t172 + t163 * t168;
	t132 = -t159 * t171 - t161 * t168;
	t1 = [0, t140 * t173 + t152 * t167, 0, t134 * t153, 0, -t135 * t152 - t140 * t153; 0, t137 * t173 + t152 * t168, 0, t132 * t153, 0, -t133 * t152 - t137 * t153; 0, -t145 * t173 + t146 * t152, 0, t142 * t153, 0, -t143 * t152 + t145 * t153; 0, -t140 * t174 + t153 * t167, 0, -t134 * t152, 0, -t135 * t153 + t140 * t152; 0, -t137 * t174 + t153 * t168, 0, -t132 * t152, 0, -t133 * t153 + t137 * t152; 0, t145 * t174 + t146 * t153, 0, -t142 * t152, 0, -t143 * t153 - t145 * t152; 0, t140 * t161, 0, t135, 0, 0; 0, t137 * t161, 0, t133, 0, 0; 0, -t145 * t161, 0, t143, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end