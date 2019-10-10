% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:53
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPPRR11_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR11_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_jacobig_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t55 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t55, 0, 0, 0, 0; 0, -cos(qJ(1)) * t55, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:51
	% EndTime: 2019-10-10 09:53:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t61 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t61, 0, 0, 0, 0; 0, -cos(qJ(1)) * t61, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:51
	% EndTime: 2019-10-10 09:53:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->3), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t75 = cos(pkin(6));
	t76 = sin(qJ(2));
	t80 = t75 * t76;
	t79 = cos(qJ(1));
	t78 = cos(qJ(2));
	t77 = sin(qJ(1));
	t74 = sin(pkin(6));
	t1 = [0, t77 * t74, 0, 0, -t77 * t80 + t79 * t78, 0; 0, -t79 * t74, 0, 0, t77 * t78 + t79 * t80, 0; 1, t75, 0, 0, t74 * t76, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:51
	% EndTime: 2019-10-10 09:53:51
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (14->9), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->16)
	t109 = sin(pkin(6));
	t112 = sin(qJ(1));
	t120 = t112 * t109;
	t111 = sin(qJ(2));
	t119 = t112 * t111;
	t113 = cos(qJ(2));
	t118 = t112 * t113;
	t114 = cos(qJ(1));
	t117 = t114 * t109;
	t116 = t114 * t111;
	t115 = t114 * t113;
	t110 = cos(pkin(6));
	t108 = pkin(11) + qJ(5);
	t107 = cos(t108);
	t106 = sin(t108);
	t1 = [0, t120, 0, 0, -t110 * t119 + t115, t106 * t120 - (t110 * t118 + t116) * t107; 0, -t117, 0, 0, t110 * t116 + t118, -t106 * t117 - (-t110 * t115 + t119) * t107; 1, t110, 0, 0, t109 * t111, t109 * t113 * t107 + t110 * t106;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end