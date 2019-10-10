% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR8
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPRR8_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR8_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobig_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(12)) * t18, 0, 0, 0, 0; 0, -cos(pkin(12)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
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
	% StartTime: 2019-10-09 22:39:10
	% EndTime: 2019-10-09 22:39:10
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->6), mult. (17->14), div. (0->0), fcn. (28->8), ass. (0->12)
	t101 = sin(pkin(6));
	t99 = sin(pkin(12));
	t109 = t99 * t101;
	t102 = cos(pkin(12));
	t108 = t102 * t101;
	t104 = cos(pkin(6));
	t106 = cos(qJ(2));
	t107 = t104 * t106;
	t105 = sin(qJ(2));
	t103 = cos(pkin(7));
	t100 = sin(pkin(7));
	t1 = [0, t109, -(-t102 * t105 - t99 * t107) * t100 + t103 * t109, 0, 0, 0; 0, -t108, -(t102 * t107 - t99 * t105) * t100 - t103 * t108, 0, 0, 0; 0, t104, -t101 * t106 * t100 + t104 * t103, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:10
	% EndTime: 2019-10-09 22:39:10
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (16->14), mult. (48->33), div. (0->0), fcn. (72->10), ass. (0->18)
	t120 = sin(pkin(12));
	t122 = sin(pkin(6));
	t134 = t120 * t122;
	t121 = sin(pkin(7));
	t133 = t121 * t122;
	t123 = cos(pkin(12));
	t132 = t123 * t122;
	t125 = cos(pkin(6));
	t127 = sin(qJ(2));
	t131 = t125 * t127;
	t129 = cos(qJ(2));
	t130 = t125 * t129;
	t128 = cos(qJ(3));
	t126 = sin(qJ(3));
	t124 = cos(pkin(7));
	t119 = -t120 * t130 - t123 * t127;
	t118 = -t120 * t127 + t123 * t130;
	t1 = [0, t134, -t119 * t121 + t124 * t134, 0, (-t120 * t131 + t123 * t129) * t128 + (t119 * t124 + t120 * t133) * t126, 0; 0, -t132, -t118 * t121 - t124 * t132, 0, (t120 * t129 + t123 * t131) * t128 + (t118 * t124 - t121 * t132) * t126, 0; 0, t125, t125 * t124 - t129 * t133, 0, t125 * t121 * t126 + (t124 * t126 * t129 + t127 * t128) * t122, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:10
	% EndTime: 2019-10-09 22:39:10
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (34->21), mult. (100->47), div. (0->0), fcn. (145->12), ass. (0->29)
	t162 = sin(pkin(12));
	t164 = sin(pkin(6));
	t182 = t162 * t164;
	t163 = sin(pkin(7));
	t181 = t163 * t164;
	t167 = cos(pkin(6));
	t180 = t163 * t167;
	t165 = cos(pkin(12));
	t179 = t165 * t164;
	t166 = cos(pkin(7));
	t173 = cos(qJ(2));
	t178 = t166 * t173;
	t170 = sin(qJ(2));
	t177 = t167 * t170;
	t176 = t167 * t173;
	t158 = -t162 * t170 + t165 * t176;
	t175 = -t158 * t166 + t163 * t179;
	t160 = -t162 * t176 - t165 * t170;
	t174 = t160 * t166 + t162 * t181;
	t172 = cos(qJ(3));
	t171 = cos(qJ(5));
	t169 = sin(qJ(3));
	t168 = sin(qJ(5));
	t161 = -t162 * t177 + t165 * t173;
	t159 = t162 * t173 + t165 * t177;
	t157 = t167 * t166 - t173 * t181;
	t156 = -t160 * t163 + t166 * t182;
	t155 = -t158 * t163 - t166 * t179;
	t1 = [0, t182, t156, 0, t161 * t172 + t174 * t169, t156 * t168 - (t161 * t169 - t174 * t172) * t171; 0, -t179, t155, 0, t159 * t172 - t175 * t169, t155 * t168 - (t159 * t169 + t175 * t172) * t171; 0, t167, t157, 0, t169 * t180 + (t169 * t178 + t170 * t172) * t164, t157 * t168 - (-t172 * t180 + (t169 * t170 - t172 * t178) * t164) * t171;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end