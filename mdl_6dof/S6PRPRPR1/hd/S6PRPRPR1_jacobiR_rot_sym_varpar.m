% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:30
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRPR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t22 = cos(pkin(6));
	t23 = sin(qJ(2));
	t26 = t22 * t23;
	t24 = cos(qJ(2));
	t25 = t22 * t24;
	t21 = cos(pkin(10));
	t20 = sin(pkin(6));
	t19 = sin(pkin(10));
	t1 = [0, -t19 * t25 - t21 * t23, 0, 0, 0, 0; 0, -t19 * t23 + t21 * t25, 0, 0, 0, 0; 0, t20 * t24, 0, 0, 0, 0; 0, t19 * t26 - t21 * t24, 0, 0, 0, 0; 0, -t19 * t24 - t21 * t26, 0, 0, 0, 0; 0, -t20 * t23, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (14->7), mult. (40->16), div. (0->0), fcn. (60->8), ass. (0->13)
	t35 = sin(pkin(11));
	t38 = cos(pkin(11));
	t41 = sin(qJ(2));
	t42 = cos(qJ(2));
	t34 = -t42 * t35 - t41 * t38;
	t33 = t41 * t35 - t42 * t38;
	t40 = cos(pkin(6));
	t39 = cos(pkin(10));
	t37 = sin(pkin(6));
	t36 = sin(pkin(10));
	t32 = t34 * t40;
	t31 = t33 * t40;
	t1 = [0, t36 * t31 + t39 * t34, 0, 0, 0, 0; 0, -t39 * t31 + t36 * t34, 0, 0, 0, 0; 0, -t33 * t37, 0, 0, 0, 0; 0, -t36 * t32 + t39 * t33, 0, 0, 0, 0; 0, t39 * t32 + t36 * t33, 0, 0, 0, 0; 0, t34 * t37, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (44->15), mult. (122->36), div. (0->0), fcn. (178->10), ass. (0->23)
	t93 = sin(pkin(6));
	t97 = sin(qJ(4));
	t104 = t93 * t97;
	t99 = cos(qJ(4));
	t103 = t93 * t99;
	t100 = cos(qJ(2));
	t91 = sin(pkin(11));
	t94 = cos(pkin(11));
	t98 = sin(qJ(2));
	t102 = t100 * t94 - t98 * t91;
	t101 = t100 * t91 + t98 * t94;
	t96 = cos(pkin(6));
	t87 = t101 * t96;
	t92 = sin(pkin(10));
	t95 = cos(pkin(10));
	t81 = t102 * t92 + t95 * t87;
	t83 = t102 * t95 - t92 * t87;
	t86 = t102 * t96;
	t85 = t101 * t93;
	t84 = t102 * t93;
	t82 = -t101 * t95 - t92 * t86;
	t80 = -t101 * t92 + t95 * t86;
	t1 = [0, t82 * t99, 0, t92 * t103 - t83 * t97, 0, 0; 0, t80 * t99, 0, -t95 * t103 - t81 * t97, 0, 0; 0, t84 * t99, 0, -t85 * t97 + t96 * t99, 0, 0; 0, -t82 * t97, 0, -t92 * t104 - t83 * t99, 0, 0; 0, -t80 * t97, 0, t95 * t104 - t81 * t99, 0, 0; 0, -t84 * t97, 0, -t85 * t99 - t96 * t97, 0, 0; 0, t83, 0, 0, 0, 0; 0, t81, 0, 0, 0, 0; 0, t85, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (62->16), mult. (122->36), div. (0->0), fcn. (178->10), ass. (0->24)
	t104 = sin(pkin(10));
	t105 = sin(pkin(6));
	t113 = t104 * t105;
	t107 = cos(pkin(10));
	t112 = t105 * t107;
	t108 = cos(pkin(6));
	t103 = sin(pkin(11));
	t106 = cos(pkin(11));
	t109 = sin(qJ(2));
	t110 = cos(qJ(2));
	t111 = t110 * t103 + t109 * t106;
	t96 = t111 * t108;
	t97 = t109 * t103 - t110 * t106;
	t90 = -t104 * t97 + t107 * t96;
	t92 = -t104 * t96 - t107 * t97;
	t102 = qJ(4) + pkin(12);
	t101 = cos(t102);
	t100 = sin(t102);
	t95 = t97 * t108;
	t94 = t111 * t105;
	t93 = t97 * t105;
	t91 = t104 * t95 - t107 * t111;
	t89 = -t104 * t111 - t107 * t95;
	t1 = [0, t91 * t101, 0, -t92 * t100 + t101 * t113, 0, 0; 0, t89 * t101, 0, -t90 * t100 - t101 * t112, 0, 0; 0, -t93 * t101, 0, -t94 * t100 + t108 * t101, 0, 0; 0, -t91 * t100, 0, -t100 * t113 - t92 * t101, 0, 0; 0, -t89 * t100, 0, t100 * t112 - t90 * t101, 0, 0; 0, t93 * t100, 0, -t108 * t100 - t94 * t101, 0, 0; 0, t92, 0, 0, 0, 0; 0, t90, 0, 0, 0, 0; 0, t94, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (153->29), mult. (319->65), div. (0->0), fcn. (454->12), ass. (0->34)
	t156 = sin(pkin(11));
	t159 = cos(pkin(11));
	t163 = sin(qJ(2));
	t165 = cos(qJ(2));
	t149 = t163 * t156 - t165 * t159;
	t155 = qJ(4) + pkin(12);
	t154 = cos(t155);
	t162 = sin(qJ(6));
	t175 = t154 * t162;
	t164 = cos(qJ(6));
	t174 = t154 * t164;
	t157 = sin(pkin(10));
	t158 = sin(pkin(6));
	t173 = t157 * t158;
	t160 = cos(pkin(10));
	t172 = t158 * t160;
	t161 = cos(pkin(6));
	t167 = t165 * t156 + t163 * t159;
	t148 = t167 * t161;
	t169 = t160 * t148 - t157 * t149;
	t168 = -t157 * t148 - t160 * t149;
	t166 = t149 * t161;
	t153 = sin(t155);
	t147 = t167 * t158;
	t146 = t149 * t158;
	t144 = t147 * t154 + t161 * t153;
	t143 = -t147 * t153 + t161 * t154;
	t141 = t157 * t166 - t160 * t167;
	t138 = -t157 * t167 - t160 * t166;
	t136 = t153 * t173 + t154 * t168;
	t135 = -t153 * t168 + t154 * t173;
	t134 = -t153 * t172 + t154 * t169;
	t133 = -t153 * t169 - t154 * t172;
	t1 = [0, t141 * t174 + t162 * t168, 0, t135 * t164, 0, -t136 * t162 - t141 * t164; 0, t138 * t174 + t162 * t169, 0, t133 * t164, 0, -t134 * t162 - t138 * t164; 0, -t146 * t174 + t147 * t162, 0, t143 * t164, 0, -t144 * t162 + t146 * t164; 0, -t141 * t175 + t164 * t168, 0, -t135 * t162, 0, -t136 * t164 + t141 * t162; 0, -t138 * t175 + t164 * t169, 0, -t133 * t162, 0, -t134 * t164 + t138 * t162; 0, t146 * t175 + t147 * t164, 0, -t143 * t162, 0, -t144 * t164 - t146 * t162; 0, t141 * t153, 0, t136, 0, 0; 0, t138 * t153, 0, t134, 0, 0; 0, -t146 * t153, 0, t144, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end