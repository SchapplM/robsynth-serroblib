% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:00
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRPR8_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR8_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_jacobig_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:04
	% EndTime: 2019-10-09 23:00:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:04
	% EndTime: 2019-10-09 23:00:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:04
	% EndTime: 2019-10-09 23:00:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(12)) * t18, 0, 0, 0, 0; 0, -cos(pkin(12)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:05
	% EndTime: 2019-10-09 23:00:05
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
	% StartTime: 2019-10-09 23:00:05
	% EndTime: 2019-10-09 23:00:05
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
	% StartTime: 2019-10-09 23:00:05
	% EndTime: 2019-10-09 23:00:05
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (16->14), mult. (48->33), div. (0->0), fcn. (72->10), ass. (0->18)
	t140 = sin(pkin(12));
	t142 = sin(pkin(6));
	t154 = t140 * t142;
	t141 = sin(pkin(7));
	t153 = t141 * t142;
	t143 = cos(pkin(12));
	t152 = t143 * t142;
	t145 = cos(pkin(6));
	t147 = sin(qJ(2));
	t151 = t145 * t147;
	t149 = cos(qJ(2));
	t150 = t145 * t149;
	t148 = cos(qJ(3));
	t146 = sin(qJ(3));
	t144 = cos(pkin(7));
	t139 = -t140 * t150 - t143 * t147;
	t138 = -t140 * t147 + t143 * t150;
	t1 = [0, t154, -t139 * t141 + t144 * t154, (-t140 * t151 + t143 * t149) * t146 + (-t139 * t144 - t140 * t153) * t148, 0, 0; 0, -t152, -t138 * t141 - t144 * t152, (t140 * t149 + t143 * t151) * t146 + (-t138 * t144 + t141 * t152) * t148, 0, 0; 0, t145, t145 * t144 - t149 * t153, -t145 * t141 * t148 + (-t144 * t148 * t149 + t146 * t147) * t142, 0, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:05
	% EndTime: 2019-10-09 23:00:05
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (34->21), mult. (100->47), div. (0->0), fcn. (145->12), ass. (0->29)
	t164 = sin(pkin(12));
	t166 = sin(pkin(6));
	t184 = t164 * t166;
	t165 = sin(pkin(7));
	t183 = t165 * t166;
	t169 = cos(pkin(6));
	t182 = t165 * t169;
	t167 = cos(pkin(12));
	t181 = t167 * t166;
	t168 = cos(pkin(7));
	t175 = cos(qJ(2));
	t180 = t168 * t175;
	t172 = sin(qJ(2));
	t179 = t169 * t172;
	t178 = t169 * t175;
	t160 = -t164 * t172 + t167 * t178;
	t177 = -t160 * t168 + t165 * t181;
	t162 = -t164 * t178 - t167 * t172;
	t176 = t162 * t168 + t164 * t183;
	t174 = cos(qJ(3));
	t173 = cos(qJ(4));
	t171 = sin(qJ(3));
	t170 = sin(qJ(4));
	t163 = -t164 * t179 + t167 * t175;
	t161 = t164 * t175 + t167 * t179;
	t159 = t169 * t168 - t175 * t183;
	t158 = -t162 * t165 + t168 * t184;
	t157 = -t160 * t165 - t168 * t181;
	t1 = [0, t184, t158, t163 * t171 - t176 * t174, 0, (t163 * t174 + t176 * t171) * t173 + t158 * t170; 0, -t181, t157, t161 * t171 + t177 * t174, 0, (t161 * t174 - t177 * t171) * t173 + t157 * t170; 0, t169, t159, -t174 * t182 + (t171 * t172 - t174 * t180) * t166, 0, (t171 * t182 + (t171 * t180 + t172 * t174) * t166) * t173 + t159 * t170;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end