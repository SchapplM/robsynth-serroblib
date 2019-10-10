% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRPPRR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:28
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPPRR3_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR3_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobig_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(10)) * t18, 0, 0, 0, 0; 0, -cos(pkin(10)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t45 = sin(pkin(6));
	t1 = [0, sin(pkin(10)) * t45, 0, 0, 0, 0; 0, -cos(pkin(10)) * t45, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t34 = sin(pkin(6));
	t1 = [0, sin(pkin(10)) * t34, 0, 0, 0, 0; 0, -cos(pkin(10)) * t34, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (8->8), mult. (22->19), div. (0->0), fcn. (35->8), ass. (0->11)
	t79 = cos(pkin(6));
	t80 = sin(qJ(2));
	t83 = t79 * t80;
	t81 = cos(qJ(2));
	t82 = t79 * t81;
	t78 = cos(pkin(10));
	t77 = cos(pkin(11));
	t76 = sin(pkin(6));
	t75 = sin(pkin(10));
	t74 = sin(pkin(11));
	t1 = [0, t75 * t76, 0, 0, (-t75 * t83 + t78 * t81) * t74 - (t75 * t82 + t78 * t80) * t77, 0; 0, -t78 * t76, 0, 0, (t75 * t81 + t78 * t83) * t74 - (t75 * t80 - t78 * t82) * t77, 0; 0, t79, 0, 0, (t74 * t80 + t77 * t81) * t76, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:09
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (18->14), mult. (50->32), div. (0->0), fcn. (76->10), ass. (0->19)
	t115 = sin(pkin(10));
	t116 = sin(pkin(6));
	t127 = t115 * t116;
	t118 = cos(pkin(10));
	t126 = t118 * t116;
	t119 = cos(pkin(6));
	t121 = sin(qJ(2));
	t125 = t119 * t121;
	t123 = cos(qJ(2));
	t124 = t119 * t123;
	t122 = cos(qJ(5));
	t120 = sin(qJ(5));
	t117 = cos(pkin(11));
	t114 = sin(pkin(11));
	t113 = -t115 * t125 + t118 * t123;
	t112 = t115 * t124 + t118 * t121;
	t111 = t115 * t123 + t118 * t125;
	t110 = t115 * t121 - t118 * t124;
	t1 = [0, t127, 0, 0, -t112 * t117 + t113 * t114, (t112 * t114 + t113 * t117) * t120 + t122 * t127; 0, -t126, 0, 0, -t110 * t117 + t111 * t114, (t110 * t114 + t111 * t117) * t120 - t122 * t126; 0, t119, 0, 0, (t114 * t121 + t117 * t123) * t116, t119 * t122 + (-t114 * t123 + t117 * t121) * t120 * t116;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end