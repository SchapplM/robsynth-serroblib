% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PPRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:18
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PPRRRR1_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRR1_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_jacobig_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (15->13), div. (0->0), fcn. (23->8), ass. (0->11)
	t41 = sin(pkin(6));
	t44 = cos(pkin(7));
	t47 = t41 * t44;
	t42 = cos(pkin(13));
	t45 = cos(pkin(6));
	t46 = t42 * t45;
	t43 = cos(pkin(12));
	t40 = sin(pkin(7));
	t39 = sin(pkin(12));
	t38 = sin(pkin(13));
	t1 = [0, 0, -(-t43 * t38 - t39 * t46) * t40 + t39 * t47, 0, 0, 0; 0, 0, -(-t39 * t38 + t43 * t46) * t40 - t43 * t47, 0, 0, 0; 0, 0, -t41 * t42 * t40 + t45 * t44, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:21
	% EndTime: 2019-10-09 21:18:21
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (15->13), mult. (46->32), div. (0->0), fcn. (67->10), ass. (0->17)
	t89 = sin(pkin(12));
	t95 = cos(pkin(6));
	t101 = t89 * t95;
	t90 = sin(pkin(7));
	t91 = sin(pkin(6));
	t100 = t90 * t91;
	t94 = cos(pkin(7));
	t99 = t91 * t94;
	t93 = cos(pkin(12));
	t98 = t93 * t95;
	t97 = cos(qJ(3));
	t96 = sin(qJ(3));
	t92 = cos(pkin(13));
	t88 = sin(pkin(13));
	t87 = -t92 * t101 - t93 * t88;
	t86 = -t89 * t88 + t92 * t98;
	t1 = [0, 0, -t87 * t90 + t89 * t99, (-t88 * t101 + t93 * t92) * t96 + (-t89 * t100 - t87 * t94) * t97, 0, 0; 0, 0, -t86 * t90 - t93 * t99, (t88 * t98 + t89 * t92) * t96 + (t93 * t100 - t86 * t94) * t97, 0, 0; 0, 0, -t92 * t100 + t95 * t94, -t95 * t90 * t97 + (-t92 * t94 * t97 + t88 * t96) * t91, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:21
	% EndTime: 2019-10-09 21:18:21
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (25->13), mult. (77->32), div. (0->0), fcn. (111->10), ass. (0->20)
	t108 = sin(pkin(12));
	t114 = cos(pkin(6));
	t120 = t108 * t114;
	t109 = sin(pkin(7));
	t110 = sin(pkin(6));
	t119 = t109 * t110;
	t113 = cos(pkin(7));
	t118 = t110 * t113;
	t112 = cos(pkin(12));
	t117 = t112 * t114;
	t116 = cos(qJ(3));
	t115 = sin(qJ(3));
	t111 = cos(pkin(13));
	t107 = sin(pkin(13));
	t106 = -t112 * t107 - t111 * t120;
	t105 = -t108 * t107 + t111 * t117;
	t104 = -t114 * t109 * t116 + (-t111 * t113 * t116 + t107 * t115) * t110;
	t103 = (-t107 * t120 + t112 * t111) * t115 + (-t106 * t113 - t108 * t119) * t116;
	t102 = (t107 * t117 + t108 * t111) * t115 + (-t105 * t113 + t112 * t119) * t116;
	t1 = [0, 0, -t106 * t109 + t108 * t118, t103, t103, 0; 0, 0, -t105 * t109 - t112 * t118, t102, t102, 0; 0, 0, -t111 * t119 + t114 * t113, t104, t104, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:21
	% EndTime: 2019-10-09 21:18:21
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (49->21), mult. (129->46), div. (0->0), fcn. (184->12), ass. (0->32)
	t166 = sin(pkin(12));
	t172 = cos(pkin(6));
	t182 = t166 * t172;
	t167 = sin(pkin(7));
	t168 = sin(pkin(6));
	t181 = t167 * t168;
	t180 = t167 * t172;
	t171 = cos(pkin(7));
	t179 = t168 * t171;
	t169 = cos(pkin(13));
	t178 = t169 * t171;
	t170 = cos(pkin(12));
	t177 = t170 * t172;
	t165 = sin(pkin(13));
	t158 = -t166 * t165 + t169 * t177;
	t176 = -t158 * t171 + t170 * t181;
	t160 = -t170 * t165 - t169 * t182;
	t175 = t160 * t171 + t166 * t181;
	t174 = cos(qJ(3));
	t173 = sin(qJ(3));
	t164 = qJ(4) + qJ(5);
	t163 = cos(t164);
	t162 = sin(t164);
	t161 = -t165 * t182 + t170 * t169;
	t159 = t165 * t177 + t166 * t169;
	t157 = -t169 * t181 + t172 * t171;
	t156 = -t160 * t167 + t166 * t179;
	t155 = -t158 * t167 - t170 * t179;
	t154 = -t174 * t180 + (t165 * t173 - t174 * t178) * t168;
	t153 = t161 * t173 - t175 * t174;
	t152 = t159 * t173 + t176 * t174;
	t1 = [0, 0, t156, t153, t153, (t161 * t174 + t175 * t173) * t162 - t156 * t163; 0, 0, t155, t152, t152, (t159 * t174 - t176 * t173) * t162 - t155 * t163; 0, 0, t157, t154, t154, (t173 * t180 + (t165 * t174 + t173 * t178) * t168) * t162 - t157 * t163;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end