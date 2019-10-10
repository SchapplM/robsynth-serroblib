% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:26
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPPRR2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
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
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
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
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (14->6), mult. (40->16), div. (0->0), fcn. (60->8), ass. (0->13)
	t59 = sin(pkin(11));
	t62 = cos(pkin(11));
	t65 = sin(qJ(2));
	t66 = cos(qJ(2));
	t67 = t66 * t59 + t65 * t62;
	t57 = t65 * t59 - t66 * t62;
	t64 = cos(pkin(6));
	t63 = cos(pkin(10));
	t61 = sin(pkin(6));
	t60 = sin(pkin(10));
	t56 = t67 * t64;
	t55 = t57 * t64;
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, -t60 * t55 + t63 * t67, 0, 0, 0, 0; 0, t63 * t55 + t60 * t67, 0, 0, 0, 0; 0, t57 * t61, 0, 0, 0, 0; 0, -t60 * t56 - t63 * t57, 0, 0, 0, 0; 0, t63 * t56 - t60 * t57, 0, 0, 0, 0; 0, t67 * t61, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (41->13), mult. (122->36), div. (0->0), fcn. (178->10), ass. (0->23)
	t95 = sin(pkin(6));
	t99 = sin(qJ(5));
	t106 = t95 * t99;
	t101 = cos(qJ(5));
	t105 = t101 * t95;
	t100 = sin(qJ(2));
	t102 = cos(qJ(2));
	t93 = sin(pkin(11));
	t96 = cos(pkin(11));
	t104 = t100 * t96 + t102 * t93;
	t91 = t100 * t93 - t102 * t96;
	t98 = cos(pkin(6));
	t103 = t91 * t98;
	t97 = cos(pkin(10));
	t94 = sin(pkin(10));
	t90 = t104 * t98;
	t89 = t104 * t95;
	t88 = t91 * t95;
	t86 = t94 * t103 - t104 * t97;
	t85 = -t94 * t90 - t97 * t91;
	t84 = -t97 * t103 - t104 * t94;
	t83 = t97 * t90 - t94 * t91;
	t1 = [0, t85 * t99, 0, 0, -t86 * t101 - t94 * t106, 0; 0, t83 * t99, 0, 0, -t84 * t101 + t97 * t106, 0; 0, t89 * t99, 0, 0, t88 * t101 - t98 * t99, 0; 0, t85 * t101, 0, 0, -t94 * t105 + t86 * t99, 0; 0, t83 * t101, 0, 0, t97 * t105 + t84 * t99, 0; 0, t89 * t101, 0, 0, -t98 * t101 - t88 * t99, 0; 0, t86, 0, 0, 0, 0; 0, t84, 0, 0, 0, 0; 0, -t88, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (117->30), mult. (319->65), div. (0->0), fcn. (454->12), ass. (0->33)
	t147 = sin(pkin(11));
	t150 = cos(pkin(11));
	t155 = sin(qJ(2));
	t158 = cos(qJ(2));
	t143 = t155 * t147 - t158 * t150;
	t149 = sin(pkin(6));
	t154 = sin(qJ(5));
	t168 = t149 * t154;
	t157 = cos(qJ(5));
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
	t138 = t140 * t154 + t152 * t157;
	t137 = t140 * t157 - t152 * t154;
	t135 = t148 * t159 - t151 * t160;
	t132 = -t148 * t160 - t151 * t159;
	t130 = -t132 * t154 - t151 * t167;
	t129 = -t132 * t157 + t151 * t168;
	t128 = -t135 * t154 + t148 * t167;
	t127 = -t135 * t157 - t148 * t168;
	t1 = [0, t135 * t153 + t161 * t165, 0, 0, t127 * t156, -t128 * t153 + t156 * t161; 0, t132 * t153 + t162 * t165, 0, 0, t129 * t156, -t130 * t153 + t156 * t162; 0, -t140 * t153 + t141 * t165, 0, 0, t137 * t156, -t138 * t153 + t141 * t156; 0, t135 * t156 - t161 * t166, 0, 0, -t127 * t153, -t128 * t156 - t153 * t161; 0, t132 * t156 - t162 * t166, 0, 0, -t129 * t153, -t130 * t156 - t153 * t162; 0, -t140 * t156 - t141 * t166, 0, 0, -t137 * t153, -t138 * t156 - t141 * t153; 0, -t161 * t157, 0, 0, t128, 0; 0, -t162 * t157, 0, 0, t130, 0; 0, -t141 * t157, 0, 0, t138, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end