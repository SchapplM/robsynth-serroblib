% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR7
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:58
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRPR7_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR7_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_jacobig_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(12)) * t18, 0, 0, 0, 0; 0, -cos(pkin(12)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (17->14), div. (0->0), fcn. (28->8), ass. (0->12)
	t64 = sin(pkin(12));
	t66 = sin(pkin(6));
	t74 = t64 * t66;
	t67 = cos(pkin(12));
	t73 = t67 * t66;
	t69 = cos(pkin(6));
	t71 = cos(qJ(2));
	t72 = t69 * t71;
	t70 = sin(qJ(2));
	t68 = cos(pkin(7));
	t65 = sin(pkin(7));
	t1 = [0, t74, -(-t64 * t72 - t67 * t70) * t65 + t68 * t74, 0, 0, 0; 0, -t73, -(-t64 * t70 + t67 * t72) * t65 - t68 * t73, 0, 0, 0; 0, t69, -t66 * t71 * t65 + t69 * t68, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (16->14), mult. (48->33), div. (0->0), fcn. (72->10), ass. (0->18)
	t119 = sin(pkin(12));
	t121 = sin(pkin(6));
	t133 = t119 * t121;
	t120 = sin(pkin(7));
	t132 = t120 * t121;
	t122 = cos(pkin(12));
	t131 = t122 * t121;
	t124 = cos(pkin(6));
	t126 = sin(qJ(2));
	t130 = t124 * t126;
	t128 = cos(qJ(2));
	t129 = t124 * t128;
	t127 = cos(qJ(3));
	t125 = sin(qJ(3));
	t123 = cos(pkin(7));
	t118 = -t119 * t129 - t122 * t126;
	t117 = -t119 * t126 + t122 * t129;
	t1 = [0, t133, -t118 * t120 + t123 * t133, (-t119 * t130 + t122 * t128) * t125 + (-t118 * t123 - t119 * t132) * t127, 0, 0; 0, -t131, -t117 * t120 - t123 * t131, (t119 * t128 + t122 * t130) * t125 + (-t117 * t123 + t120 * t131) * t127, 0, 0; 0, t124, t124 * t123 - t128 * t132, -t124 * t120 * t127 + (-t123 * t127 * t128 + t125 * t126) * t121, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (16->14), mult. (48->33), div. (0->0), fcn. (72->10), ass. (0->18)
	t154 = sin(pkin(12));
	t156 = sin(pkin(6));
	t168 = t154 * t156;
	t155 = sin(pkin(7));
	t167 = t155 * t156;
	t157 = cos(pkin(12));
	t166 = t157 * t156;
	t159 = cos(pkin(6));
	t161 = sin(qJ(2));
	t165 = t159 * t161;
	t163 = cos(qJ(2));
	t164 = t159 * t163;
	t162 = cos(qJ(3));
	t160 = sin(qJ(3));
	t158 = cos(pkin(7));
	t153 = -t154 * t164 - t157 * t161;
	t152 = -t154 * t161 + t157 * t164;
	t1 = [0, t168, -t153 * t155 + t158 * t168, (-t154 * t165 + t157 * t163) * t160 + (-t153 * t158 - t154 * t167) * t162, 0, 0; 0, -t166, -t152 * t155 - t158 * t166, (t154 * t163 + t157 * t165) * t160 + (-t152 * t158 + t155 * t166) * t162, 0, 0; 0, t159, t159 * t158 - t163 * t167, -t159 * t155 * t162 + (-t158 * t162 * t163 + t160 * t161) * t156, 0, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:05
	% EndTime: 2019-10-09 22:58:05
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (34->21), mult. (100->47), div. (0->0), fcn. (145->12), ass. (0->29)
	t170 = sin(pkin(12));
	t172 = sin(pkin(6));
	t190 = t170 * t172;
	t171 = sin(pkin(7));
	t189 = t171 * t172;
	t175 = cos(pkin(6));
	t188 = t171 * t175;
	t173 = cos(pkin(12));
	t187 = t173 * t172;
	t174 = cos(pkin(7));
	t181 = cos(qJ(2));
	t186 = t174 * t181;
	t178 = sin(qJ(2));
	t185 = t175 * t178;
	t184 = t175 * t181;
	t166 = -t170 * t178 + t173 * t184;
	t183 = -t166 * t174 + t171 * t187;
	t168 = -t170 * t184 - t173 * t178;
	t182 = t168 * t174 + t170 * t189;
	t180 = cos(qJ(3));
	t179 = cos(qJ(4));
	t177 = sin(qJ(3));
	t176 = sin(qJ(4));
	t169 = -t170 * t185 + t173 * t181;
	t167 = t170 * t181 + t173 * t185;
	t165 = t175 * t174 - t181 * t189;
	t164 = -t168 * t171 + t174 * t190;
	t163 = -t166 * t171 - t174 * t187;
	t1 = [0, t190, t164, t169 * t177 - t182 * t180, 0, (t169 * t180 + t182 * t177) * t176 - t164 * t179; 0, -t187, t163, t167 * t177 + t183 * t180, 0, (t167 * t180 - t183 * t177) * t176 - t163 * t179; 0, t175, t165, -t180 * t188 + (t177 * t178 - t180 * t186) * t172, 0, (t177 * t188 + (t177 * t186 + t178 * t180) * t172) * t176 - t165 * t179;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end