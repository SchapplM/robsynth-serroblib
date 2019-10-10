% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR1
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:23
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRRPR1_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR1_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->1), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->2)
	t20 = qJ(1) + pkin(10);
	t1 = [0, 0, qJD(1) * cos(t20), 0, 0, 0; 0, 0, qJD(1) * sin(t20), 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t34 = qJ(1) + pkin(10);
	t33 = qJD(1) * cos(t34);
	t32 = qJD(1) * sin(t34);
	t1 = [0, 0, t33, t33, 0, 0; 0, 0, t32, t32, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t39 = qJ(1) + pkin(10);
	t38 = qJD(1) * cos(t39);
	t37 = qJD(1) * sin(t39);
	t1 = [0, 0, t38, t38, 0, 0; 0, 0, t37, t37, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:32
	% EndTime: 2019-10-10 01:23:32
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (23->5), mult. (13->8), div. (0->0), fcn. (13->4), ass. (0->10)
	t128 = qJ(3) + qJ(4) + pkin(11);
	t129 = qJD(3) + qJD(4);
	t131 = cos(t128) * t129;
	t130 = qJ(1) + pkin(10);
	t126 = sin(t130);
	t122 = qJD(1) * t126;
	t127 = cos(t130);
	t123 = qJD(1) * t127;
	t124 = sin(t128);
	t1 = [0, 0, t123, t123, 0, -t124 * t122 + t127 * t131; 0, 0, t122, t122, 0, t124 * t123 + t126 * t131; 0, 0, 0, 0, 0, t129 * t124;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end