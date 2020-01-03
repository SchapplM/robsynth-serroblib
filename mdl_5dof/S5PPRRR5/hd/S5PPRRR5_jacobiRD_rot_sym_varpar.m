% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PPRRR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:35:59
	% EndTime: 2019-12-31 17:35:59
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:35:59
	% EndTime: 2019-12-31 17:35:59
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:35:59
	% EndTime: 2019-12-31 17:35:59
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:35:59
	% EndTime: 2019-12-31 17:35:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->3), mult. (16->6), div. (0->0), fcn. (16->4), ass. (0->7)
	t55 = cos(qJ(3));
	t54 = sin(qJ(3));
	t53 = cos(pkin(8));
	t52 = sin(pkin(8));
	t51 = (t52 * t54 + t53 * t55) * qJD(3);
	t50 = (-t52 * t55 + t53 * t54) * qJD(3);
	t1 = [0, 0, -t51, 0, 0; 0, 0, t50, 0, 0; 0, 0, 0, 0, 0; 0, 0, t50, 0, 0; 0, 0, t51, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:35:59
	% EndTime: 2019-12-31 17:35:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (40->6), mult. (32->6), div. (0->0), fcn. (32->4), ass. (0->9)
	t83 = cos(pkin(8));
	t82 = sin(pkin(8));
	t81 = qJ(3) + qJ(4);
	t80 = qJD(3) + qJD(4);
	t79 = cos(t81);
	t78 = sin(t81);
	t77 = (t78 * t82 + t79 * t83) * t80;
	t76 = (t78 * t83 - t79 * t82) * t80;
	t1 = [0, 0, -t77, -t77, 0; 0, 0, t76, t76, 0; 0, 0, 0, 0, 0; 0, 0, t76, t76, 0; 0, 0, t77, t77, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:35:59
	% EndTime: 2019-12-31 17:35:59
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (128->15), mult. (126->20), div. (0->0), fcn. (138->6), ass. (0->21)
	t118 = sin(pkin(8));
	t117 = qJ(3) + qJ(4);
	t108 = sin(qJ(5));
	t116 = qJD(5) * t108;
	t109 = cos(qJ(5));
	t115 = qJD(5) * t109;
	t114 = qJD(3) + qJD(4);
	t113 = cos(t117);
	t106 = sin(t117);
	t112 = t106 * t114;
	t107 = cos(pkin(8));
	t101 = -t118 * t106 - t107 * t113;
	t110 = t114 * t113;
	t99 = t107 * t112 - t118 * t110;
	t111 = t101 * t115 + t99 * t108;
	t96 = -t101 * t116 + t99 * t109;
	t102 = -t107 * t106 + t118 * t113;
	t100 = -t107 * t110 - t118 * t112;
	t98 = t100 * t109 - t102 * t116;
	t97 = -t100 * t108 - t102 * t115;
	t1 = [0, 0, t98, t98, t111; 0, 0, t96, t96, t97; 0, 0, 0, 0, t116; 0, 0, t97, t97, t96; 0, 0, -t111, -t111, -t98; 0, 0, 0, 0, t115; 0, 0, -t99, -t99, 0; 0, 0, t100, t100, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end