% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRPRR4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:52:16
	% EndTime: 2019-12-05 15:52:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:52:16
	% EndTime: 2019-12-05 15:52:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:52:16
	% EndTime: 2019-12-05 15:52:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t22 = cos(pkin(5));
	t23 = sin(qJ(2));
	t26 = t22 * t23;
	t24 = cos(qJ(2));
	t25 = t22 * t24;
	t21 = cos(pkin(9));
	t20 = sin(pkin(5));
	t19 = sin(pkin(9));
	t1 = [0, -t19 * t25 - t21 * t23, 0, 0, 0; 0, -t19 * t23 + t21 * t25, 0, 0, 0; 0, t20 * t24, 0, 0, 0; 0, t19 * t26 - t21 * t24, 0, 0, 0; 0, -t19 * t24 - t21 * t26, 0, 0, 0; 0, -t20 * t23, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:52:16
	% EndTime: 2019-12-05 15:52:16
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (14->7), mult. (40->16), div. (0->0), fcn. (60->8), ass. (0->13)
	t35 = sin(pkin(10));
	t38 = cos(pkin(10));
	t41 = sin(qJ(2));
	t42 = cos(qJ(2));
	t34 = -t42 * t35 - t41 * t38;
	t33 = t41 * t35 - t42 * t38;
	t40 = cos(pkin(5));
	t39 = cos(pkin(9));
	t37 = sin(pkin(5));
	t36 = sin(pkin(9));
	t32 = t34 * t40;
	t31 = t33 * t40;
	t1 = [0, t36 * t31 + t39 * t34, 0, 0, 0; 0, -t39 * t31 + t36 * t34, 0, 0, 0; 0, -t33 * t37, 0, 0, 0; 0, -t36 * t32 + t39 * t33, 0, 0, 0; 0, t39 * t32 + t36 * t33, 0, 0, 0; 0, t34 * t37, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:52:16
	% EndTime: 2019-12-05 15:52:16
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (44->15), mult. (122->36), div. (0->0), fcn. (178->10), ass. (0->23)
	t93 = sin(pkin(5));
	t97 = sin(qJ(4));
	t104 = t93 * t97;
	t99 = cos(qJ(4));
	t103 = t93 * t99;
	t100 = cos(qJ(2));
	t91 = sin(pkin(10));
	t94 = cos(pkin(10));
	t98 = sin(qJ(2));
	t102 = t100 * t94 - t98 * t91;
	t101 = t100 * t91 + t98 * t94;
	t96 = cos(pkin(5));
	t87 = t101 * t96;
	t92 = sin(pkin(9));
	t95 = cos(pkin(9));
	t81 = t102 * t92 + t95 * t87;
	t83 = t102 * t95 - t92 * t87;
	t86 = t102 * t96;
	t85 = t101 * t93;
	t84 = t102 * t93;
	t82 = -t101 * t95 - t92 * t86;
	t80 = -t101 * t92 + t95 * t86;
	t1 = [0, t82 * t99, 0, t92 * t103 - t83 * t97, 0; 0, t80 * t99, 0, -t95 * t103 - t81 * t97, 0; 0, t84 * t99, 0, -t85 * t97 + t96 * t99, 0; 0, -t82 * t97, 0, -t92 * t104 - t83 * t99, 0; 0, -t80 * t97, 0, t95 * t104 - t81 * t99, 0; 0, -t84 * t97, 0, -t85 * t99 - t96 * t97, 0; 0, t83, 0, 0, 0; 0, t81, 0, 0, 0; 0, t85, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:52:17
	% EndTime: 2019-12-05 15:52:17
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (114->28), mult. (319->65), div. (0->0), fcn. (454->12), ass. (0->33)
	t144 = sin(pkin(10));
	t147 = cos(pkin(10));
	t152 = sin(qJ(2));
	t155 = cos(qJ(2));
	t140 = t152 * t144 - t155 * t147;
	t146 = sin(pkin(5));
	t151 = sin(qJ(4));
	t165 = t146 * t151;
	t154 = cos(qJ(4));
	t164 = t146 * t154;
	t150 = sin(qJ(5));
	t163 = t150 * t154;
	t153 = cos(qJ(5));
	t161 = t153 * t154;
	t149 = cos(pkin(5));
	t157 = t155 * t144 + t152 * t147;
	t139 = t157 * t149;
	t145 = sin(pkin(9));
	t148 = cos(pkin(9));
	t159 = t148 * t139 - t145 * t140;
	t158 = -t145 * t139 - t148 * t140;
	t156 = t140 * t149;
	t138 = t157 * t146;
	t137 = t140 * t146;
	t135 = t138 * t154 + t149 * t151;
	t134 = -t138 * t151 + t149 * t154;
	t132 = t145 * t156 - t148 * t157;
	t129 = -t145 * t157 - t148 * t156;
	t127 = t145 * t165 + t154 * t158;
	t126 = t145 * t164 - t151 * t158;
	t125 = -t148 * t165 + t154 * t159;
	t124 = -t148 * t164 - t151 * t159;
	t1 = [0, t132 * t161 + t150 * t158, 0, t126 * t153, -t127 * t150 - t132 * t153; 0, t129 * t161 + t150 * t159, 0, t124 * t153, -t125 * t150 - t129 * t153; 0, -t137 * t161 + t138 * t150, 0, t134 * t153, -t135 * t150 + t137 * t153; 0, -t132 * t163 + t153 * t158, 0, -t126 * t150, -t127 * t153 + t132 * t150; 0, -t129 * t163 + t153 * t159, 0, -t124 * t150, -t125 * t153 + t129 * t150; 0, t137 * t163 + t138 * t153, 0, -t134 * t150, -t135 * t153 - t137 * t150; 0, t132 * t151, 0, t127, 0; 0, t129 * t151, 0, t125, 0; 0, -t137 * t151, 0, t135, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end