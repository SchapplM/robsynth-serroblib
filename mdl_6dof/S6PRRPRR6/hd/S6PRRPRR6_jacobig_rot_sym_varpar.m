% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:35
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPRR6_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR6_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobig_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(12)) * t18, 0, 0, 0, 0; 0, -cos(pkin(12)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
	% DurationCPUTime: 0.03s
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
	% StartTime: 2019-10-09 22:35:20
	% EndTime: 2019-10-09 22:35:20
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->6), mult. (17->14), div. (0->0), fcn. (28->8), ass. (0->12)
	t117 = sin(pkin(12));
	t119 = sin(pkin(6));
	t127 = t117 * t119;
	t120 = cos(pkin(12));
	t126 = t120 * t119;
	t122 = cos(pkin(6));
	t124 = cos(qJ(2));
	t125 = t122 * t124;
	t123 = sin(qJ(2));
	t121 = cos(pkin(7));
	t118 = sin(pkin(7));
	t1 = [0, t127, -(-t117 * t125 - t120 * t123) * t118 + t121 * t127, 0, 0, 0; 0, -t126, -(-t117 * t123 + t120 * t125) * t118 - t121 * t126, 0, 0, 0; 0, t122, -t119 * t124 * t118 + t122 * t121, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:20
	% EndTime: 2019-10-09 22:35:20
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (16->14), mult. (48->33), div. (0->0), fcn. (72->10), ass. (0->18)
	t126 = sin(pkin(12));
	t128 = sin(pkin(6));
	t140 = t126 * t128;
	t127 = sin(pkin(7));
	t139 = t127 * t128;
	t129 = cos(pkin(12));
	t138 = t129 * t128;
	t131 = cos(pkin(6));
	t133 = sin(qJ(2));
	t137 = t131 * t133;
	t135 = cos(qJ(2));
	t136 = t131 * t135;
	t134 = cos(qJ(3));
	t132 = sin(qJ(3));
	t130 = cos(pkin(7));
	t125 = -t126 * t136 - t129 * t133;
	t124 = -t126 * t133 + t129 * t136;
	t1 = [0, t140, -t125 * t127 + t130 * t140, 0, (-t126 * t137 + t129 * t135) * t132 + (-t125 * t130 - t126 * t139) * t134, 0; 0, -t138, -t124 * t127 - t130 * t138, 0, (t126 * t135 + t129 * t137) * t132 + (-t124 * t130 + t127 * t138) * t134, 0; 0, t131, t131 * t130 - t135 * t139, 0, -t131 * t127 * t134 + (-t130 * t134 * t135 + t132 * t133) * t128, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:20
	% EndTime: 2019-10-09 22:35:20
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (40->22), mult. (100->47), div. (0->0), fcn. (145->12), ass. (0->30)
	t171 = sin(pkin(12));
	t173 = sin(pkin(6));
	t189 = t171 * t173;
	t172 = sin(pkin(7));
	t188 = t172 * t173;
	t176 = cos(pkin(6));
	t187 = t172 * t176;
	t174 = cos(pkin(12));
	t186 = t174 * t173;
	t175 = cos(pkin(7));
	t180 = cos(qJ(2));
	t185 = t175 * t180;
	t178 = sin(qJ(2));
	t184 = t176 * t178;
	t183 = t176 * t180;
	t164 = -t171 * t178 + t174 * t183;
	t182 = -t164 * t175 + t172 * t186;
	t166 = -t171 * t183 - t174 * t178;
	t181 = t166 * t175 + t171 * t188;
	t179 = cos(qJ(3));
	t177 = sin(qJ(3));
	t170 = pkin(13) + qJ(5);
	t169 = cos(t170);
	t168 = sin(t170);
	t167 = -t171 * t184 + t174 * t180;
	t165 = t171 * t180 + t174 * t184;
	t163 = t176 * t175 - t180 * t188;
	t162 = -t166 * t172 + t175 * t189;
	t161 = -t164 * t172 - t175 * t186;
	t1 = [0, t189, t162, 0, t167 * t177 - t181 * t179, (t167 * t179 + t181 * t177) * t168 - t162 * t169; 0, -t186, t161, 0, t165 * t177 + t182 * t179, (t165 * t179 - t182 * t177) * t168 - t161 * t169; 0, t176, t163, 0, -t179 * t187 + (t177 * t178 - t179 * t185) * t173, (t177 * t187 + (t177 * t185 + t178 * t179) * t173) * t168 - t163 * t169;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end