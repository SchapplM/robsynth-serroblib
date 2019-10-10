% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:04
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPRPRR11_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR11_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_jacobig_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (5->5), mult. (15->13), div. (0->0), fcn. (23->8), ass. (0->11)
	t78 = sin(pkin(6));
	t80 = cos(pkin(7));
	t85 = t78 * t80;
	t79 = cos(pkin(12));
	t81 = cos(pkin(6));
	t84 = t79 * t81;
	t83 = cos(qJ(1));
	t82 = sin(qJ(1));
	t77 = sin(pkin(7));
	t76 = sin(pkin(12));
	t1 = [0, 0, -(-t83 * t76 - t82 * t84) * t77 + t82 * t85, 0, 0, 0; 0, 0, -(-t82 * t76 + t83 * t84) * t77 - t83 * t85, 0, 0, 0; 1, 0, -t78 * t79 * t77 + t81 * t80, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (15->13), div. (0->0), fcn. (23->8), ass. (0->11)
	t114 = sin(pkin(6));
	t116 = cos(pkin(7));
	t121 = t114 * t116;
	t115 = cos(pkin(12));
	t117 = cos(pkin(6));
	t120 = t115 * t117;
	t119 = cos(qJ(1));
	t118 = sin(qJ(1));
	t113 = sin(pkin(7));
	t112 = sin(pkin(12));
	t1 = [0, 0, -(-t119 * t112 - t118 * t120) * t113 + t118 * t121, 0, 0, 0; 0, 0, -(-t118 * t112 + t119 * t120) * t113 - t119 * t121, 0, 0, 0; 1, 0, -t114 * t115 * t113 + t117 * t116, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (15->13), mult. (46->31), div. (0->0), fcn. (67->10), ass. (0->19)
	t128 = sin(pkin(6));
	t133 = sin(qJ(1));
	t141 = t128 * t133;
	t135 = cos(qJ(1));
	t140 = t128 * t135;
	t126 = sin(pkin(12));
	t139 = t133 * t126;
	t129 = cos(pkin(12));
	t138 = t133 * t129;
	t137 = t135 * t126;
	t136 = t135 * t129;
	t134 = cos(qJ(3));
	t132 = sin(qJ(3));
	t131 = cos(pkin(6));
	t130 = cos(pkin(7));
	t127 = sin(pkin(7));
	t125 = -t131 * t138 - t137;
	t124 = t131 * t136 - t139;
	t1 = [0, 0, -t125 * t127 + t130 * t141, 0, (-t131 * t139 + t136) * t132 + (-t125 * t130 - t127 * t141) * t134, 0; 0, 0, -t124 * t127 - t130 * t140, 0, (t131 * t137 + t138) * t132 + (-t124 * t130 + t127 * t140) * t134, 0; 1, 0, -t128 * t129 * t127 + t131 * t130, 0, -t131 * t127 * t134 + (-t129 * t130 * t134 + t126 * t132) * t128, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:01
	% EndTime: 2019-10-10 01:04:01
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (39->21), mult. (98->45), div. (0->0), fcn. (140->12), ass. (0->31)
	t174 = sin(pkin(7));
	t178 = cos(pkin(6));
	t192 = t174 * t178;
	t175 = sin(pkin(6));
	t180 = sin(qJ(1));
	t191 = t175 * t180;
	t182 = cos(qJ(1));
	t190 = t175 * t182;
	t176 = cos(pkin(12));
	t177 = cos(pkin(7));
	t189 = t176 * t177;
	t173 = sin(pkin(12));
	t188 = t180 * t173;
	t187 = t180 * t176;
	t186 = t182 * t173;
	t185 = t182 * t176;
	t166 = t178 * t185 - t188;
	t184 = -t166 * t177 + t174 * t190;
	t168 = -t178 * t187 - t186;
	t183 = t168 * t177 + t174 * t191;
	t181 = cos(qJ(3));
	t179 = sin(qJ(3));
	t172 = pkin(13) + qJ(5);
	t171 = cos(t172);
	t170 = sin(t172);
	t169 = -t178 * t188 + t185;
	t167 = t178 * t186 + t187;
	t165 = -t175 * t176 * t174 + t178 * t177;
	t164 = -t168 * t174 + t177 * t191;
	t163 = -t166 * t174 - t177 * t190;
	t1 = [0, 0, t164, 0, t169 * t179 - t183 * t181, (t169 * t181 + t183 * t179) * t170 - t164 * t171; 0, 0, t163, 0, t167 * t179 + t184 * t181, (t167 * t181 - t184 * t179) * t170 - t163 * t171; 1, 0, t165, 0, -t181 * t192 + (t173 * t179 - t181 * t189) * t175, (t179 * t192 + (t173 * t181 + t179 * t189) * t175) * t170 - t165 * t171;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end