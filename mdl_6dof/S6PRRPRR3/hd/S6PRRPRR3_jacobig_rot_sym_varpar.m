% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR3
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
% Datum: 2019-10-09 22:29
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPRR3_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR3_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobig_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(12)) * t18, 0, 0, 0, 0; 0, -cos(pkin(12)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
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
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (17->14), div. (0->0), fcn. (28->8), ass. (0->12)
	t85 = sin(pkin(12));
	t87 = sin(pkin(6));
	t95 = t85 * t87;
	t88 = cos(pkin(12));
	t94 = t88 * t87;
	t90 = cos(pkin(6));
	t92 = cos(qJ(2));
	t93 = t90 * t92;
	t91 = sin(qJ(2));
	t89 = cos(pkin(7));
	t86 = sin(pkin(7));
	t1 = [0, t95, -(-t85 * t93 - t88 * t91) * t86 + t89 * t95, 0, 0, 0; 0, -t94, -(-t85 * t91 + t88 * t93) * t86 - t89 * t94, 0, 0, 0; 0, t90, -t87 * t92 * t86 + t90 * t89, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (25->16), mult. (72->35), div. (0->0), fcn. (105->12), ass. (0->23)
	t133 = sin(pkin(12));
	t135 = sin(pkin(6));
	t148 = t133 * t135;
	t137 = cos(pkin(12));
	t147 = t137 * t135;
	t139 = cos(pkin(6));
	t141 = sin(qJ(2));
	t146 = t139 * t141;
	t143 = cos(qJ(2));
	t145 = t139 * t143;
	t132 = sin(pkin(13));
	t136 = cos(pkin(13));
	t140 = sin(qJ(3));
	t142 = cos(qJ(3));
	t144 = -t132 * t140 + t136 * t142;
	t138 = cos(pkin(7));
	t134 = sin(pkin(7));
	t131 = -t142 * t132 - t140 * t136;
	t130 = -t133 * t145 - t137 * t141;
	t129 = -t133 * t141 + t137 * t145;
	t128 = t144 * t138;
	t127 = t144 * t134;
	t1 = [0, t148, -t130 * t134 + t138 * t148, 0, -(-t133 * t146 + t137 * t143) * t131 - t130 * t128 - t127 * t148, 0; 0, -t147, -t129 * t134 - t138 * t147, 0, -(t133 * t143 + t137 * t146) * t131 - t129 * t128 + t127 * t147, 0; 0, t139, -t135 * t143 * t134 + t139 * t138, 0, -t139 * t127 + (-t128 * t143 - t131 * t141) * t135, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:27
	% EndTime: 2019-10-09 22:29:27
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (52->25), mult. (148->53), div. (0->0), fcn. (211->14), ass. (0->32)
	t180 = sin(pkin(12));
	t182 = sin(pkin(6));
	t197 = t180 * t182;
	t184 = cos(pkin(12));
	t196 = t184 * t182;
	t186 = cos(pkin(6));
	t189 = sin(qJ(2));
	t195 = t186 * t189;
	t192 = cos(qJ(2));
	t194 = t186 * t192;
	t179 = sin(pkin(13));
	t183 = cos(pkin(13));
	t188 = sin(qJ(3));
	t191 = cos(qJ(3));
	t193 = t191 * t179 + t188 * t183;
	t178 = -t188 * t179 + t191 * t183;
	t190 = cos(qJ(5));
	t187 = sin(qJ(5));
	t185 = cos(pkin(7));
	t181 = sin(pkin(7));
	t176 = -t180 * t195 + t184 * t192;
	t175 = -t180 * t194 - t184 * t189;
	t174 = t180 * t192 + t184 * t195;
	t173 = -t180 * t189 + t184 * t194;
	t172 = -t182 * t192 * t181 + t186 * t185;
	t171 = t193 * t185;
	t170 = t178 * t185;
	t169 = t193 * t181;
	t168 = t178 * t181;
	t167 = -t175 * t181 + t185 * t197;
	t166 = -t173 * t181 - t185 * t196;
	t1 = [0, t197, t167, 0, -t168 * t197 - t175 * t170 + t176 * t193, (t169 * t197 + t175 * t171 + t176 * t178) * t187 - t167 * t190; 0, -t196, t166, 0, t168 * t196 - t173 * t170 + t174 * t193, (-t169 * t196 + t173 * t171 + t174 * t178) * t187 - t166 * t190; 0, t186, t172, 0, -t186 * t168 + (-t170 * t192 + t189 * t193) * t182, (t186 * t169 + (t171 * t192 + t178 * t189) * t182) * t187 - t172 * t190;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end