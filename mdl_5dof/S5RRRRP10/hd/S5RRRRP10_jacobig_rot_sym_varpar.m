% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% Jg_rot [3x5]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S5RRRRP10_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_jacobig_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP10_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_jacobig_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:44:20
	% EndTime: 2019-12-29 20:44:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:44:20
	% EndTime: 2019-12-29 20:44:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:44:20
	% EndTime: 2019-12-29 20:44:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(5));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0; 1, cos(pkin(5)), 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:44:20
	% EndTime: 2019-12-29 20:44:20
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t70 = cos(pkin(5));
	t73 = cos(qJ(2));
	t75 = t70 * t73;
	t74 = cos(qJ(1));
	t72 = sin(qJ(1));
	t71 = sin(qJ(2));
	t69 = sin(pkin(5));
	t1 = [0, t72 * t69, t74 * t71 + t72 * t75, 0, 0; 0, -t74 * t69, t72 * t71 - t74 * t75, 0, 0; 1, t70, -t69 * t73, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:44:21
	% EndTime: 2019-12-29 20:44:21
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (9->9), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->15)
	t101 = sin(pkin(5));
	t105 = sin(qJ(1));
	t114 = t105 * t101;
	t104 = sin(qJ(2));
	t113 = t105 * t104;
	t107 = cos(qJ(2));
	t112 = t105 * t107;
	t108 = cos(qJ(1));
	t111 = t108 * t101;
	t110 = t108 * t104;
	t109 = t108 * t107;
	t106 = cos(qJ(3));
	t103 = sin(qJ(3));
	t102 = cos(pkin(5));
	t1 = [0, t114, t102 * t112 + t110, (-t102 * t113 + t109) * t103 - t106 * t114, 0; 0, -t111, -t102 * t109 + t113, (t102 * t110 + t112) * t103 + t106 * t111, 0; 1, t102, -t101 * t107, t101 * t104 * t103 - t102 * t106, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:44:20
	% EndTime: 2019-12-29 20:44:21
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (9->9), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->15)
	t103 = sin(pkin(5));
	t107 = sin(qJ(1));
	t116 = t107 * t103;
	t106 = sin(qJ(2));
	t115 = t107 * t106;
	t109 = cos(qJ(2));
	t114 = t107 * t109;
	t110 = cos(qJ(1));
	t113 = t110 * t103;
	t112 = t110 * t106;
	t111 = t110 * t109;
	t108 = cos(qJ(3));
	t105 = sin(qJ(3));
	t104 = cos(pkin(5));
	t1 = [0, t116, t104 * t114 + t112, (-t104 * t115 + t111) * t105 - t108 * t116, 0; 0, -t113, -t104 * t111 + t115, (t104 * t112 + t114) * t105 + t108 * t113, 0; 1, t104, -t103 * t109, t103 * t106 * t105 - t104 * t108, 0;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,5);
end