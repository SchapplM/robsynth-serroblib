% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRP6
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% JgD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S5RRRRP6_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_jacobigD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_jacobigD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP6_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_jacobigD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:33:16
	% EndTime: 2019-12-29 20:33:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:33:16
	% EndTime: 2019-12-29 20:33:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:33:09
	% EndTime: 2019-12-29 20:33:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:33:14
	% EndTime: 2019-12-29 20:33:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t30 = qJD(1) * cos(qJ(1));
	t29 = qJD(1) * sin(qJ(1));
	t1 = [0, t30, t30, 0, 0; 0, t29, t29, 0, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:33:10
	% EndTime: 2019-12-29 20:33:10
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (10->4), mult. (13->8), div. (0->0), fcn. (13->4), ass. (0->9)
	t116 = qJD(2) + qJD(3);
	t117 = qJ(2) + qJ(3);
	t120 = cos(t117) * t116;
	t118 = sin(qJ(1));
	t112 = qJD(1) * t118;
	t119 = cos(qJ(1));
	t113 = qJD(1) * t119;
	t114 = sin(t117);
	t1 = [0, t113, t113, -t114 * t112 + t119 * t120, 0; 0, t112, t112, t114 * t113 + t118 * t120, 0; 0, 0, 0, t116 * t114, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:33:14
	% EndTime: 2019-12-29 20:33:14
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (10->4), mult. (13->8), div. (0->0), fcn. (13->4), ass. (0->9)
	t114 = qJD(2) + qJD(3);
	t115 = qJ(2) + qJ(3);
	t118 = cos(t115) * t114;
	t116 = sin(qJ(1));
	t110 = qJD(1) * t116;
	t117 = cos(qJ(1));
	t111 = qJD(1) * t117;
	t112 = sin(t115);
	t1 = [0, t111, t111, -t112 * t110 + t117 * t118, 0; 0, t110, t110, t112 * t111 + t116 * t118, 0; 0, 0, 0, t114 * t112, 0;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,5);
end