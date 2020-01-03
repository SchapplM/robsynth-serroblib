% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% JRD_rot [9x4]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S4RPPR6_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_jacobiRD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_jacobiRD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPPR6_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_jacobiRD_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:40:54
	% EndTime: 2019-12-31 16:40:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:40:54
	% EndTime: 2019-12-31 16:40:54
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0; -t31, 0, 0, 0; 0, 0, 0, 0; t31, 0, 0, 0; -t30, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:40:54
	% EndTime: 2019-12-31 16:40:54
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t22 = qJD(1) * sin(qJ(1));
	t21 = qJD(1) * cos(qJ(1));
	t18 = cos(pkin(6));
	t17 = sin(pkin(6));
	t1 = [-t18 * t21, 0, 0, 0; -t18 * t22, 0, 0, 0; 0, 0, 0, 0; t17 * t21, 0, 0, 0; t17 * t22, 0, 0, 0; 0, 0, 0, 0; -t22, 0, 0, 0; t21, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:40:54
	% EndTime: 2019-12-31 16:40:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->5), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t104 = qJD(1) * sin(qJ(1));
	t103 = qJD(1) * cos(qJ(1));
	t100 = cos(pkin(6));
	t99 = sin(pkin(6));
	t1 = [-t100 * t103, 0, 0, 0; -t100 * t104, 0, 0, 0; 0, 0, 0, 0; -t104, 0, 0, 0; t103, 0, 0, 0; 0, 0, 0, 0; -t99 * t103, 0, 0, 0; -t99 * t104, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:40:54
	% EndTime: 2019-12-31 16:40:54
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (27->9), mult. (106->18), div. (0->0), fcn. (106->6), ass. (0->18)
	t57 = sin(qJ(1));
	t66 = qJD(1) * t57;
	t54 = sin(pkin(6));
	t55 = cos(pkin(6));
	t56 = sin(qJ(4));
	t58 = cos(qJ(4));
	t65 = -t54 * t58 + t55 * t56;
	t64 = t54 * t56 + t55 * t58;
	t59 = cos(qJ(1));
	t63 = t64 * t59;
	t62 = qJD(1) * t65;
	t61 = t65 * qJD(4);
	t60 = t64 * qJD(4);
	t53 = -qJD(1) * t63 + t57 * t61;
	t52 = t57 * t60 + t59 * t62;
	t51 = t59 * t61 + t64 * t66;
	t50 = -qJD(4) * t63 + t57 * t62;
	t1 = [t53, 0, 0, t50; -t51, 0, 0, -t52; 0, 0, 0, t61; t52, 0, 0, t51; t50, 0, 0, t53; 0, 0, 0, t60; t66, 0, 0, 0; -qJD(1) * t59, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,4);
end