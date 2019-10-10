% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPPP1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:15
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPP1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPP1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:12
	% EndTime: 2019-10-10 11:15:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:12
	% EndTime: 2019-10-10 11:15:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:12
	% EndTime: 2019-10-10 11:15:12
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t10 = sin(qJ(1));
	t9 = sin(qJ(2));
	t16 = t10 * t9;
	t12 = cos(qJ(1));
	t15 = t12 * t9;
	t11 = cos(qJ(2));
	t14 = t10 * t11;
	t13 = t12 * t11;
	t1 = [-t14, -t15, 0, 0, 0, 0; t13, -t16, 0, 0, 0, 0; 0, t11, 0, 0, 0, 0; t16, -t13, 0, 0, 0, 0; -t15, -t14, 0, 0, 0, 0; 0, -t9, 0, 0, 0, 0; t12, 0, 0, 0, 0, 0; t10, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:12
	% EndTime: 2019-10-10 11:15:12
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (14->12), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
	t60 = sin(qJ(2));
	t61 = sin(qJ(1));
	t71 = t61 * t60;
	t62 = cos(qJ(3));
	t70 = t61 * t62;
	t59 = sin(qJ(3));
	t63 = cos(qJ(2));
	t69 = t63 * t59;
	t68 = t63 * t62;
	t64 = cos(qJ(1));
	t67 = t64 * t60;
	t66 = t64 * t62;
	t65 = t64 * t63;
	t58 = t61 * t59 + t62 * t65;
	t57 = -t59 * t65 + t70;
	t56 = t64 * t59 - t61 * t68;
	t55 = t61 * t69 + t66;
	t1 = [t56, -t60 * t66, t57, 0, 0, 0; t58, -t60 * t70, -t55, 0, 0, 0; 0, t68, -t60 * t59, 0, 0, 0; t55, t59 * t67, -t58, 0, 0, 0; t57, t59 * t71, t56, 0, 0, 0; 0, -t69, -t60 * t62, 0, 0, 0; -t71, t65, 0, 0, 0, 0; t67, t61 * t63, 0, 0, 0, 0; 0, t60, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:12
	% EndTime: 2019-10-10 11:15:13
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (53->27), mult. (183->71), div. (0->0), fcn. (258->10), ass. (0->35)
	t100 = sin(pkin(10));
	t103 = cos(pkin(6));
	t129 = t100 * t103;
	t101 = sin(pkin(6));
	t105 = sin(qJ(2));
	t128 = t101 * t105;
	t102 = cos(pkin(10));
	t127 = t102 * t103;
	t107 = cos(qJ(3));
	t126 = t103 * t107;
	t104 = sin(qJ(3));
	t108 = cos(qJ(2));
	t125 = t104 * t108;
	t124 = t105 * t103;
	t106 = sin(qJ(1));
	t123 = t105 * t106;
	t122 = t105 * t107;
	t109 = cos(qJ(1));
	t121 = t105 * t109;
	t120 = t106 * t104;
	t119 = t107 * t108;
	t118 = t109 * t104;
	t117 = t109 * t107;
	t96 = t108 * t120 + t117;
	t116 = -t101 * t123 + t103 * t96;
	t98 = t106 * t107 - t108 * t118;
	t115 = t101 * t121 + t103 * t98;
	t114 = -t103 * t125 + t128;
	t113 = t101 * t108 + t104 * t124;
	t112 = t103 * t108 - t104 * t128;
	t111 = t100 * t122 + t113 * t102;
	t110 = t113 * t100 - t102 * t122;
	t99 = t108 * t117 + t120;
	t97 = -t106 * t119 + t118;
	t1 = [t116 * t100 + t97 * t102, t110 * t109, t98 * t102 - t99 * t129, 0, 0, 0; t115 * t100 + t99 * t102, t110 * t106, -t96 * t102 + t97 * t129, 0, 0, 0; 0, t114 * t100 + t102 * t119, (-t100 * t126 - t102 * t104) * t105, 0, 0, 0; -t97 * t100 + t116 * t102, t111 * t109, -t98 * t100 - t99 * t127, 0, 0, 0; -t99 * t100 + t115 * t102, t111 * t106, t96 * t100 + t97 * t127, 0, 0, 0; 0, -t100 * t119 + t114 * t102, (t100 * t104 - t102 * t126) * t105, 0, 0, 0; -t96 * t101 - t103 * t123, t112 * t109, t99 * t101, 0, 0, 0; -t98 * t101 + t103 * t121, t112 * t106, -t97 * t101, 0, 0, 0; 0, t101 * t125 + t124, t101 * t122, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:13
	% EndTime: 2019-10-10 11:15:13
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (53->27), mult. (183->71), div. (0->0), fcn. (258->10), ass. (0->35)
	t145 = sin(pkin(10));
	t148 = cos(pkin(6));
	t174 = t145 * t148;
	t146 = sin(pkin(6));
	t150 = sin(qJ(2));
	t173 = t146 * t150;
	t147 = cos(pkin(10));
	t172 = t147 * t148;
	t152 = cos(qJ(3));
	t171 = t148 * t152;
	t149 = sin(qJ(3));
	t153 = cos(qJ(2));
	t170 = t149 * t153;
	t169 = t150 * t148;
	t151 = sin(qJ(1));
	t168 = t150 * t151;
	t167 = t150 * t152;
	t154 = cos(qJ(1));
	t166 = t150 * t154;
	t165 = t151 * t149;
	t164 = t152 * t153;
	t163 = t154 * t149;
	t162 = t154 * t152;
	t141 = t153 * t165 + t162;
	t161 = -t141 * t148 + t146 * t168;
	t143 = t151 * t152 - t153 * t163;
	t160 = -t143 * t148 - t146 * t166;
	t159 = t148 * t170 - t173;
	t158 = -t146 * t153 - t149 * t169;
	t157 = t148 * t153 - t149 * t173;
	t156 = -t145 * t167 + t158 * t147;
	t155 = t158 * t145 + t147 * t167;
	t144 = t153 * t162 + t165;
	t142 = -t151 * t164 + t163;
	t1 = [-t141 * t146 - t148 * t168, t157 * t154, t144 * t146, 0, 0, 0; -t143 * t146 + t148 * t166, t157 * t151, -t142 * t146, 0, 0, 0; 0, t146 * t170 + t169, t146 * t167, 0, 0, 0; -t142 * t147 + t161 * t145, t155 * t154, -t143 * t147 + t144 * t174, 0, 0, 0; -t144 * t147 + t160 * t145, t155 * t151, t141 * t147 - t142 * t174, 0, 0, 0; 0, t159 * t145 - t147 * t164, (t145 * t171 + t147 * t149) * t150, 0, 0, 0; t142 * t145 + t161 * t147, t156 * t154, t143 * t145 + t144 * t172, 0, 0, 0; t144 * t145 + t160 * t147, t156 * t151, -t141 * t145 - t142 * t172, 0, 0, 0; 0, t145 * t164 + t159 * t147, (-t145 * t149 + t147 * t171) * t150, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:13
	% EndTime: 2019-10-10 11:15:13
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (53->27), mult. (183->71), div. (0->0), fcn. (258->10), ass. (0->35)
	t159 = sin(pkin(10));
	t162 = cos(pkin(6));
	t188 = t159 * t162;
	t160 = sin(pkin(6));
	t164 = sin(qJ(2));
	t187 = t160 * t164;
	t161 = cos(pkin(10));
	t186 = t161 * t162;
	t166 = cos(qJ(3));
	t185 = t162 * t166;
	t163 = sin(qJ(3));
	t167 = cos(qJ(2));
	t184 = t163 * t167;
	t183 = t164 * t162;
	t165 = sin(qJ(1));
	t182 = t164 * t165;
	t181 = t164 * t166;
	t168 = cos(qJ(1));
	t180 = t164 * t168;
	t179 = t165 * t163;
	t178 = t166 * t167;
	t177 = t168 * t163;
	t176 = t168 * t166;
	t155 = t167 * t179 + t176;
	t175 = -t155 * t162 + t160 * t182;
	t157 = t165 * t166 - t167 * t177;
	t174 = t157 * t162 + t160 * t180;
	t173 = t162 * t184 - t187;
	t172 = t160 * t167 + t163 * t183;
	t171 = t162 * t167 - t163 * t187;
	t170 = -t159 * t181 - t172 * t161;
	t169 = t172 * t159 - t161 * t181;
	t158 = t167 * t176 + t179;
	t156 = -t165 * t178 + t177;
	t1 = [-t155 * t160 - t162 * t182, t171 * t168, t158 * t160, 0, 0, 0; -t157 * t160 + t162 * t180, t171 * t165, -t156 * t160, 0, 0, 0; 0, t160 * t184 + t183, t160 * t181, 0, 0, 0; t156 * t159 + t175 * t161, t170 * t168, t157 * t159 + t158 * t186, 0, 0, 0; t158 * t159 - t174 * t161, t170 * t165, -t155 * t159 - t156 * t186, 0, 0, 0; 0, t159 * t178 + t173 * t161, (-t159 * t163 + t161 * t185) * t164, 0, 0, 0; t156 * t161 - t175 * t159, t169 * t168, t157 * t161 - t158 * t188, 0, 0, 0; t158 * t161 + t174 * t159, t169 * t165, -t155 * t161 + t156 * t188, 0, 0, 0; 0, -t173 * t159 + t161 * t178, (-t159 * t185 - t161 * t163) * t164, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end