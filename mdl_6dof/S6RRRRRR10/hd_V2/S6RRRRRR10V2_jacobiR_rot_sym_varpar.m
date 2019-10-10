% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRRR10V2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:38
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR10V2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR10V2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_jacobiR_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
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
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t8 = sin(qJ(2));
	t9 = sin(qJ(1));
	t15 = t9 * t8;
	t11 = cos(qJ(1));
	t14 = t11 * t8;
	t10 = cos(qJ(2));
	t13 = t9 * t10;
	t12 = t11 * t10;
	t1 = [-t13, -t14, 0, 0, 0, 0; t12, -t15, 0, 0, 0, 0; 0, t10, 0, 0, 0, 0; t15, -t12, 0, 0, 0, 0; -t14, -t13, 0, 0, 0, 0; 0, -t8, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (28->13), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t20 = qJ(2) + qJ(3);
	t18 = sin(t20);
	t21 = sin(qJ(1));
	t26 = t21 * t18;
	t19 = cos(t20);
	t25 = t21 * t19;
	t22 = cos(qJ(1));
	t24 = t22 * t18;
	t23 = t22 * t19;
	t1 = [-t25, -t24, -t24, 0, 0, 0; t23, -t26, -t26, 0, 0, 0; 0, t19, t19, 0, 0, 0; t26, -t23, -t23, 0, 0, 0; -t24, -t25, -t25, 0, 0, 0; 0, -t18, -t18, 0, 0, 0; t22, 0, 0, 0, 0, 0; t21, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (47->16), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
	t100 = sin(qJ(4));
	t99 = qJ(2) + qJ(3);
	t98 = cos(t99);
	t110 = t98 * t100;
	t101 = sin(qJ(1));
	t109 = t101 * t100;
	t102 = cos(qJ(4));
	t108 = t101 * t102;
	t103 = cos(qJ(1));
	t107 = t103 * t100;
	t106 = t103 * t102;
	t97 = sin(t99);
	t105 = t97 * t108;
	t104 = t97 * t106;
	t96 = t103 * t98;
	t95 = t98 * t102;
	t94 = t101 * t98;
	t93 = t97 * t107;
	t92 = t97 * t109;
	t91 = t98 * t106 + t109;
	t90 = -t98 * t107 + t108;
	t89 = -t98 * t108 + t107;
	t88 = t98 * t109 + t106;
	t1 = [t89, -t104, -t104, t90, 0, 0; t91, -t105, -t105, -t88, 0, 0; 0, t95, t95, -t97 * t100, 0, 0; t88, t93, t93, -t91, 0, 0; t90, t92, t92, t89, 0, 0; 0, -t110, -t110, -t97 * t102, 0, 0; -t101 * t97, t96, t96, 0, 0, 0; t103 * t97, t94, t94, 0, 0, 0; 0, t97, t97, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (106->23), mult. (149->41), div. (0->0), fcn. (226->8), ass. (0->38)
	t131 = qJ(2) + qJ(3);
	t129 = sin(t131);
	t132 = sin(qJ(5));
	t152 = t129 * t132;
	t135 = cos(qJ(5));
	t151 = t129 * t135;
	t137 = cos(qJ(1));
	t150 = t129 * t137;
	t136 = cos(qJ(4));
	t149 = t132 * t136;
	t133 = sin(qJ(4));
	t134 = sin(qJ(1));
	t148 = t134 * t133;
	t147 = t134 * t136;
	t146 = t135 * t136;
	t145 = t137 * t133;
	t144 = t137 * t136;
	t143 = t129 * t148;
	t142 = t129 * t145;
	t130 = cos(t131);
	t123 = t130 * t147 - t145;
	t141 = -t123 * t132 + t134 * t151;
	t140 = -t123 * t135 - t134 * t152;
	t139 = -t129 * t146 + t130 * t132;
	t138 = t129 * t149 + t130 * t135;
	t126 = t130 * t133;
	t125 = t130 * t144 + t148;
	t124 = t130 * t145 - t147;
	t122 = -t130 * t148 - t144;
	t121 = t130 * t146 + t152;
	t120 = -t130 * t149 + t151;
	t119 = t139 * t137;
	t118 = t138 * t137;
	t117 = t139 * t134;
	t116 = t138 * t134;
	t115 = t125 * t135 + t132 * t150;
	t114 = -t125 * t132 + t135 * t150;
	t1 = [t140, t119, t119, -t124 * t135, t114, 0; t115, t117, t117, t122 * t135, t141, 0; 0, t121, t121, -t133 * t151, -t138, 0; -t141, t118, t118, t124 * t132, -t115, 0; t114, t116, t116, -t122 * t132, t140, 0; 0, t120, t120, t133 * t152, t139, 0; t122, -t142, -t142, t125, 0, 0; t124, -t143, -t143, t123, 0, 0; 0, t126, t126, t129 * t136, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:03
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (220->39), mult. (341->83), div. (0->0), fcn. (490->10), ass. (0->54)
	t198 = qJ(2) + qJ(3);
	t197 = cos(t198);
	t201 = sin(qJ(4));
	t206 = cos(qJ(1));
	t211 = t206 * t201;
	t202 = sin(qJ(1));
	t205 = cos(qJ(4));
	t214 = t202 * t205;
	t189 = t197 * t214 - t211;
	t204 = cos(qJ(5));
	t196 = sin(t198);
	t200 = sin(qJ(5));
	t222 = t196 * t200;
	t176 = t189 * t204 + t202 * t222;
	t210 = t206 * t205;
	t215 = t202 * t201;
	t188 = t197 * t215 + t210;
	t199 = sin(qJ(6));
	t203 = cos(qJ(6));
	t226 = t176 * t199 - t188 * t203;
	t225 = -t176 * t203 - t188 * t199;
	t221 = t196 * t204;
	t220 = t196 * t206;
	t219 = t199 * t201;
	t218 = t199 * t204;
	t217 = t200 * t205;
	t216 = t201 * t203;
	t213 = t203 * t204;
	t212 = t204 * t205;
	t209 = t196 * t219;
	t208 = t196 * t216;
	t207 = t196 * t211;
	t175 = -t189 * t200 + t202 * t221;
	t185 = t196 * t212 - t197 * t200;
	t184 = -t196 * t217 - t197 * t204;
	t191 = t197 * t210 + t215;
	t190 = t197 * t211 - t214;
	t187 = t197 * t212 + t222;
	t186 = t197 * t217 - t221;
	t183 = t185 * t206;
	t182 = t184 * t206;
	t181 = t185 * t202;
	t180 = t184 * t202;
	t179 = t191 * t204 + t200 * t220;
	t178 = t191 * t200 - t204 * t220;
	t174 = t187 * t203 + t197 * t219;
	t173 = -t187 * t199 + t197 * t216;
	t172 = -t183 * t203 - t199 * t207;
	t171 = t183 * t199 - t203 * t207;
	t170 = -t181 * t203 - t202 * t209;
	t169 = t181 * t199 - t202 * t208;
	t168 = t179 * t203 + t190 * t199;
	t167 = -t179 * t199 + t190 * t203;
	t1 = [t225, t172, t172, -t190 * t213 + t191 * t199, -t178 * t203, t167; t168, t170, t170, -t188 * t213 + t189 * t199, t175 * t203, -t226; 0, t174, t174, (t199 * t205 - t201 * t213) * t196, t184 * t203, -t185 * t199 + t208; t226, t171, t171, t190 * t218 + t191 * t203, t178 * t199, -t168; t167, t169, t169, t188 * t218 + t189 * t203, -t175 * t199, t225; 0, t173, t173, (t201 * t218 + t203 * t205) * t196, -t184 * t199, -t185 * t203 - t209; t175, t182, t182, -t190 * t200, t179, 0; t178, t180, t180, -t188 * t200, t176, 0; 0, t186, t186, -t201 * t222, t185, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end