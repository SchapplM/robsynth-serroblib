% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPRRPR9_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR9_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_jacobig_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:32
	% EndTime: 2019-10-10 01:37:32
	% DurationCPUTime: 0.08s
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
	% StartTime: 2019-10-10 01:37:32
	% EndTime: 2019-10-10 01:37:32
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (15->13), mult. (46->31), div. (0->0), fcn. (67->10), ass. (0->19)
	t121 = sin(pkin(6));
	t126 = sin(qJ(1));
	t134 = t121 * t126;
	t128 = cos(qJ(1));
	t133 = t121 * t128;
	t119 = sin(pkin(12));
	t132 = t126 * t119;
	t122 = cos(pkin(12));
	t131 = t126 * t122;
	t130 = t128 * t119;
	t129 = t128 * t122;
	t127 = cos(qJ(3));
	t125 = sin(qJ(3));
	t124 = cos(pkin(6));
	t123 = cos(pkin(7));
	t120 = sin(pkin(7));
	t118 = -t124 * t131 - t130;
	t117 = t124 * t129 - t132;
	t1 = [0, 0, -t118 * t120 + t123 * t134, (-t124 * t132 + t129) * t125 + (-t118 * t123 - t120 * t134) * t127, 0, 0; 0, 0, -t117 * t120 - t123 * t133, (t124 * t130 + t131) * t125 + (-t117 * t123 + t120 * t133) * t127, 0, 0; 1, 0, -t121 * t122 * t120 + t124 * t123, -t124 * t120 * t127 + (-t122 * t123 * t127 + t119 * t125) * t121, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:32
	% EndTime: 2019-10-10 01:37:32
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (15->13), mult. (46->31), div. (0->0), fcn. (67->10), ass. (0->19)
	t131 = sin(pkin(6));
	t136 = sin(qJ(1));
	t144 = t131 * t136;
	t138 = cos(qJ(1));
	t143 = t131 * t138;
	t129 = sin(pkin(12));
	t142 = t136 * t129;
	t132 = cos(pkin(12));
	t141 = t136 * t132;
	t140 = t138 * t129;
	t139 = t138 * t132;
	t137 = cos(qJ(3));
	t135 = sin(qJ(3));
	t134 = cos(pkin(6));
	t133 = cos(pkin(7));
	t130 = sin(pkin(7));
	t128 = -t134 * t141 - t140;
	t127 = t134 * t139 - t142;
	t1 = [0, 0, -t128 * t130 + t133 * t144, (-t134 * t142 + t139) * t135 + (-t128 * t133 - t130 * t144) * t137, 0, 0; 0, 0, -t127 * t130 - t133 * t143, (t134 * t140 + t141) * t135 + (-t127 * t133 + t130 * t143) * t137, 0, 0; 1, 0, -t131 * t132 * t130 + t134 * t133, -t134 * t130 * t137 + (-t132 * t133 * t137 + t129 * t135) * t131, 0, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:32
	% EndTime: 2019-10-10 01:37:32
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (39->21), mult. (98->45), div. (0->0), fcn. (140->12), ass. (0->31)
	t173 = sin(pkin(7));
	t177 = cos(pkin(6));
	t191 = t173 * t177;
	t174 = sin(pkin(6));
	t179 = sin(qJ(1));
	t190 = t174 * t179;
	t181 = cos(qJ(1));
	t189 = t174 * t181;
	t175 = cos(pkin(12));
	t176 = cos(pkin(7));
	t188 = t175 * t176;
	t172 = sin(pkin(12));
	t187 = t179 * t172;
	t186 = t179 * t175;
	t185 = t181 * t172;
	t184 = t181 * t175;
	t165 = t177 * t184 - t187;
	t183 = -t165 * t176 + t173 * t189;
	t167 = -t177 * t186 - t185;
	t182 = t167 * t176 + t173 * t190;
	t180 = cos(qJ(3));
	t178 = sin(qJ(3));
	t171 = qJ(4) + pkin(13);
	t170 = cos(t171);
	t169 = sin(t171);
	t168 = -t177 * t187 + t184;
	t166 = t177 * t185 + t186;
	t164 = -t174 * t175 * t173 + t177 * t176;
	t163 = -t167 * t173 + t176 * t190;
	t162 = -t165 * t173 - t176 * t189;
	t1 = [0, 0, t163, t168 * t178 - t182 * t180, 0, (t168 * t180 + t182 * t178) * t169 - t163 * t170; 0, 0, t162, t166 * t178 + t183 * t180, 0, (t166 * t180 - t183 * t178) * t169 - t162 * t170; 1, 0, t164, -t180 * t191 + (t172 * t178 - t180 * t188) * t174, 0, (t178 * t191 + (t172 * t180 + t178 * t188) * t174) * t169 - t164 * t170;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end