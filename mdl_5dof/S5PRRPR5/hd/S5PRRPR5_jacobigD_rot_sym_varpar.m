% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S5PRRPR5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% JgD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S5PRRPR5_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_jacobigD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_jacobigD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR5_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_jacobigD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:28:54
	% EndTime: 2019-12-05 16:28:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:28:54
	% EndTime: 2019-12-05 16:28:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:28:54
	% EndTime: 2019-12-05 16:28:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:28:54
	% EndTime: 2019-12-05 16:28:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (12->9), div. (0->0), fcn. (12->6), ass. (0->6)
	t82 = sin(qJ(2));
	t84 = cos(pkin(5)) * t82;
	t83 = cos(qJ(2));
	t80 = cos(pkin(9));
	t79 = sin(pkin(9));
	t1 = [0, 0, (-t79 * t84 + t80 * t83) * qJD(2), 0, 0; 0, 0, (t79 * t83 + t80 * t84) * qJD(2), 0, 0; 0, 0, sin(pkin(5)) * qJD(2) * t82, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:28:54
	% EndTime: 2019-12-05 16:28:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (12->9), div. (0->0), fcn. (12->6), ass. (0->6)
	t90 = sin(qJ(2));
	t92 = cos(pkin(5)) * t90;
	t91 = cos(qJ(2));
	t88 = cos(pkin(9));
	t87 = sin(pkin(9));
	t1 = [0, 0, (-t87 * t92 + t88 * t91) * qJD(2), 0, 0; 0, 0, (t87 * t91 + t88 * t92) * qJD(2), 0, 0; 0, 0, sin(pkin(5)) * qJD(2) * t90, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:28:55
	% EndTime: 2019-12-05 16:28:55
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (21->11), mult. (48->30), div. (0->0), fcn. (50->8), ass. (0->17)
	t139 = qJ(3) + pkin(10);
	t137 = sin(t139);
	t141 = sin(pkin(5));
	t152 = t141 * t137;
	t143 = cos(pkin(5));
	t144 = sin(qJ(2));
	t151 = t143 * t144;
	t145 = cos(qJ(2));
	t150 = t143 * t145;
	t149 = qJD(2) * t137;
	t148 = qJD(2) * t141;
	t140 = sin(pkin(9));
	t142 = cos(pkin(9));
	t147 = t140 * t145 + t142 * t151;
	t146 = -t140 * t151 + t142 * t145;
	t138 = cos(t139);
	t1 = [0, 0, t146 * qJD(2), 0, (t146 * t138 + t140 * t152) * qJD(3) + (-t140 * t150 - t142 * t144) * t149; 0, 0, t147 * qJD(2), 0, (t147 * t138 - t142 * t152) * qJD(3) + (-t140 * t144 + t142 * t150) * t149; 0, 0, t144 * t148, 0, t145 * t137 * t148 + (t138 * t141 * t144 + t137 * t143) * qJD(3);];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,5);
end