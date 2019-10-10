% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:28
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPPRR3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
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
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t49 = cos(pkin(6));
	t50 = sin(qJ(2));
	t53 = t49 * t50;
	t51 = cos(qJ(2));
	t52 = t49 * t51;
	t48 = cos(pkin(10));
	t47 = sin(pkin(6));
	t46 = sin(pkin(10));
	t1 = [0, -t46 * t52 - t48 * t50, 0, 0, 0, 0; 0, -t46 * t50 + t48 * t52, 0, 0, 0, 0; 0, t47 * t51, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, -t46 * t53 + t48 * t51, 0, 0, 0, 0; 0, t46 * t51 + t48 * t53, 0, 0, 0, 0; 0, t47 * t50, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (14->10), mult. (40->24), div. (0->0), fcn. (60->8), ass. (0->15)
	t44 = cos(pkin(6));
	t45 = sin(qJ(2));
	t48 = t44 * t45;
	t46 = cos(qJ(2));
	t47 = t44 * t46;
	t43 = cos(pkin(10));
	t42 = cos(pkin(11));
	t41 = sin(pkin(6));
	t40 = sin(pkin(10));
	t39 = sin(pkin(11));
	t38 = -t40 * t48 + t43 * t46;
	t37 = -t40 * t47 - t43 * t45;
	t36 = t40 * t46 + t43 * t48;
	t35 = -t40 * t45 + t43 * t47;
	t1 = [0, t37 * t42 + t38 * t39, 0, 0, 0, 0; 0, t35 * t42 + t36 * t39, 0, 0, 0, 0; 0, (t39 * t45 + t42 * t46) * t41, 0, 0, 0, 0; 0, -t37 * t39 + t38 * t42, 0, 0, 0, 0; 0, -t35 * t39 + t36 * t42, 0, 0, 0, 0; 0, (-t39 * t46 + t42 * t45) * t41, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (44->22), mult. (122->44), div. (0->0), fcn. (178->10), ass. (0->25)
	t101 = sin(qJ(5));
	t97 = sin(pkin(6));
	t108 = t101 * t97;
	t103 = cos(qJ(5));
	t107 = t103 * t97;
	t100 = cos(pkin(6));
	t102 = sin(qJ(2));
	t106 = t100 * t102;
	t104 = cos(qJ(2));
	t105 = t100 * t104;
	t96 = sin(pkin(10));
	t99 = cos(pkin(10));
	t90 = t96 * t102 - t99 * t105;
	t91 = t96 * t104 + t99 * t106;
	t95 = sin(pkin(11));
	t98 = cos(pkin(11));
	t85 = t90 * t95 + t91 * t98;
	t92 = t99 * t102 + t96 * t105;
	t93 = t99 * t104 - t96 * t106;
	t87 = t92 * t95 + t93 * t98;
	t89 = (t102 * t98 - t104 * t95) * t97;
	t88 = (t102 * t95 + t104 * t98) * t97;
	t86 = -t92 * t98 + t93 * t95;
	t84 = -t90 * t98 + t91 * t95;
	t1 = [0, t86 * t103, 0, 0, -t87 * t101 - t96 * t107, 0; 0, t84 * t103, 0, 0, -t85 * t101 + t99 * t107, 0; 0, t88 * t103, 0, 0, -t100 * t103 - t89 * t101, 0; 0, -t86 * t101, 0, 0, -t87 * t103 + t96 * t108, 0; 0, -t84 * t101, 0, 0, -t85 * t103 - t99 * t108, 0; 0, -t88 * t101, 0, 0, t100 * t101 - t89 * t103, 0; 0, -t87, 0, 0, 0, 0; 0, -t85, 0, 0, 0, 0; 0, -t89, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:09
	% EndTime: 2019-10-09 21:28:09
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (114->31), mult. (319->73), div. (0->0), fcn. (454->12), ass. (0->35)
	t154 = sin(pkin(6));
	t159 = sin(qJ(5));
	t171 = t154 * t159;
	t162 = cos(qJ(5));
	t170 = t154 * t162;
	t157 = cos(pkin(6));
	t160 = sin(qJ(2));
	t169 = t157 * t160;
	t163 = cos(qJ(2));
	t168 = t157 * t163;
	t158 = sin(qJ(6));
	t167 = t158 * t162;
	t161 = cos(qJ(6));
	t166 = t161 * t162;
	t153 = sin(pkin(10));
	t156 = cos(pkin(10));
	t146 = t153 * t160 - t156 * t168;
	t147 = t153 * t163 + t156 * t169;
	t152 = sin(pkin(11));
	t155 = cos(pkin(11));
	t165 = -t146 * t155 + t147 * t152;
	t148 = t153 * t168 + t156 * t160;
	t149 = -t153 * t169 + t156 * t163;
	t164 = -t148 * t155 + t149 * t152;
	t135 = t146 * t152 + t147 * t155;
	t139 = t148 * t152 + t149 * t155;
	t145 = (-t152 * t163 + t155 * t160) * t154;
	t144 = (t152 * t160 + t155 * t163) * t154;
	t141 = t145 * t162 - t157 * t159;
	t140 = -t145 * t159 - t157 * t162;
	t131 = t139 * t162 - t153 * t171;
	t130 = -t139 * t159 - t153 * t170;
	t129 = t135 * t162 + t156 * t171;
	t128 = -t135 * t159 + t156 * t170;
	t1 = [0, -t139 * t158 + t164 * t166, 0, 0, t130 * t161, -t131 * t158 + t161 * t164; 0, -t135 * t158 + t165 * t166, 0, 0, t128 * t161, -t129 * t158 + t161 * t165; 0, t144 * t166 - t145 * t158, 0, 0, t140 * t161, -t141 * t158 + t144 * t161; 0, -t139 * t161 - t164 * t167, 0, 0, -t130 * t158, -t131 * t161 - t158 * t164; 0, -t135 * t161 - t165 * t167, 0, 0, -t128 * t158, -t129 * t161 - t158 * t165; 0, -t144 * t167 - t145 * t161, 0, 0, -t140 * t158, -t141 * t161 - t144 * t158; 0, t164 * t159, 0, 0, t131, 0; 0, t165 * t159, 0, 0, t129, 0; 0, t144 * t159, 0, 0, t141, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end