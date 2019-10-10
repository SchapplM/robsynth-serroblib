% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:44
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRRP2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRP2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:44:35
	% EndTime: 2019-10-09 21:44:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:44:35
	% EndTime: 2019-10-09 21:44:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:44:35
	% EndTime: 2019-10-09 21:44:35
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
	% StartTime: 2019-10-09 21:44:35
	% EndTime: 2019-10-09 21:44:35
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
	% StartTime: 2019-10-09 21:44:36
	% EndTime: 2019-10-09 21:44:36
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
	% StartTime: 2019-10-09 21:44:36
	% EndTime: 2019-10-09 21:44:36
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (114->28), mult. (319->65), div. (0->0), fcn. (454->12), ass. (0->33)
	t144 = sin(pkin(11));
	t147 = cos(pkin(11));
	t152 = sin(qJ(2));
	t155 = cos(qJ(2));
	t140 = t152 * t144 - t155 * t147;
	t146 = sin(pkin(6));
	t151 = sin(qJ(4));
	t165 = t146 * t151;
	t154 = cos(qJ(4));
	t164 = t146 * t154;
	t150 = sin(qJ(5));
	t163 = t150 * t154;
	t153 = cos(qJ(5));
	t161 = t153 * t154;
	t149 = cos(pkin(6));
	t157 = t155 * t144 + t152 * t147;
	t139 = t157 * t149;
	t145 = sin(pkin(10));
	t148 = cos(pkin(10));
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
	t1 = [0, t132 * t161 + t150 * t158, 0, t126 * t153, -t127 * t150 - t132 * t153, 0; 0, t129 * t161 + t150 * t159, 0, t124 * t153, -t125 * t150 - t129 * t153, 0; 0, -t137 * t161 + t138 * t150, 0, t134 * t153, -t135 * t150 + t137 * t153, 0; 0, -t132 * t163 + t153 * t158, 0, -t126 * t150, -t127 * t153 + t132 * t150, 0; 0, -t129 * t163 + t153 * t159, 0, -t124 * t150, -t125 * t153 + t129 * t150, 0; 0, t137 * t163 + t138 * t153, 0, -t134 * t150, -t135 * t153 - t137 * t150, 0; 0, t132 * t151, 0, t127, 0, 0; 0, t129 * t151, 0, t125, 0, 0; 0, -t137 * t151, 0, t135, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:44:36
	% EndTime: 2019-10-09 21:44:36
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (111->25), mult. (319->65), div. (0->0), fcn. (454->12), ass. (0->33)
	t169 = sin(pkin(11));
	t172 = cos(pkin(11));
	t177 = sin(qJ(2));
	t180 = cos(qJ(2));
	t165 = t177 * t169 - t180 * t172;
	t171 = sin(pkin(6));
	t176 = sin(qJ(4));
	t190 = t171 * t176;
	t179 = cos(qJ(4));
	t189 = t171 * t179;
	t175 = sin(qJ(5));
	t188 = t175 * t179;
	t178 = cos(qJ(5));
	t186 = t178 * t179;
	t174 = cos(pkin(6));
	t182 = t180 * t169 + t177 * t172;
	t164 = t182 * t174;
	t170 = sin(pkin(10));
	t173 = cos(pkin(10));
	t184 = t173 * t164 - t170 * t165;
	t183 = -t170 * t164 - t173 * t165;
	t181 = t165 * t174;
	t163 = t182 * t171;
	t162 = t165 * t171;
	t160 = t163 * t179 + t174 * t176;
	t159 = -t163 * t176 + t174 * t179;
	t157 = t170 * t181 - t173 * t182;
	t154 = -t170 * t182 - t173 * t181;
	t152 = t170 * t190 + t179 * t183;
	t151 = t170 * t189 - t176 * t183;
	t150 = -t173 * t190 + t179 * t184;
	t149 = -t173 * t189 - t176 * t184;
	t1 = [0, t157 * t186 + t175 * t183, 0, t151 * t178, -t152 * t175 - t157 * t178, 0; 0, t154 * t186 + t175 * t184, 0, t149 * t178, -t150 * t175 - t154 * t178, 0; 0, -t162 * t186 + t163 * t175, 0, t159 * t178, -t160 * t175 + t162 * t178, 0; 0, t157 * t176, 0, t152, 0, 0; 0, t154 * t176, 0, t150, 0, 0; 0, -t162 * t176, 0, t160, 0, 0; 0, t157 * t188 - t178 * t183, 0, t151 * t175, t152 * t178 - t157 * t175, 0; 0, t154 * t188 - t178 * t184, 0, t149 * t175, t150 * t178 - t154 * t175, 0; 0, -t162 * t188 - t163 * t178, 0, t159 * t175, t160 * t178 + t162 * t175, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end