% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S4RRPR4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% JRD_rot [9x4]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 13:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S4RRPR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_jacobiRD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_jacobiRD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRPR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_jacobiRD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:33:50
	% EndTime: 2019-12-29 13:33:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:33:50
	% EndTime: 2019-12-29 13:33:50
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0; -t31, 0, 0, 0; 0, 0, 0, 0; t31, 0, 0, 0; -t30, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:33:50
	% EndTime: 2019-12-29 13:33:50
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (22->8), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t47 = qJD(1) + qJD(2);
	t48 = qJ(1) + qJ(2);
	t49 = t47 * cos(t48);
	t44 = t47 * sin(t48);
	t1 = [-t49, -t49, 0, 0; -t44, -t44, 0, 0; 0, 0, 0, 0; t44, t44, 0, 0; -t49, -t49, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:33:50
	% EndTime: 2019-12-29 13:33:50
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (30->8), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->13)
	t55 = qJ(1) + qJ(2);
	t52 = sin(t55);
	t54 = qJD(1) + qJD(2);
	t62 = t54 * t52;
	t61 = t54 * sin(pkin(7));
	t60 = t54 * cos(pkin(7));
	t59 = t52 * t60;
	t53 = cos(t55);
	t58 = t53 * t60;
	t51 = t54 * t53;
	t50 = t53 * t61;
	t49 = t52 * t61;
	t1 = [-t58, -t58, 0, 0; -t59, -t59, 0, 0; 0, 0, 0, 0; t50, t50, 0, 0; t49, t49, 0, 0; 0, 0, 0, 0; -t62, -t62, 0, 0; t51, t51, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:33:45
	% EndTime: 2019-12-29 13:33:45
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (86->14), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->16)
	t78 = qJ(1) + qJ(2);
	t74 = sin(t78);
	t77 = qJD(1) + qJD(2);
	t81 = t77 * t74;
	t75 = cos(t78);
	t71 = t77 * t75;
	t80 = qJD(4) * t74;
	t79 = qJD(4) * t75;
	t76 = pkin(7) + qJ(4);
	t73 = cos(t76);
	t72 = sin(t76);
	t70 = -t73 * t71 + t72 * t80;
	t69 = t72 * t71 + t73 * t80;
	t68 = t72 * t79 + t73 * t81;
	t67 = t72 * t81 - t73 * t79;
	t1 = [t70, t70, 0, t67; -t68, -t68, 0, -t69; 0, 0, 0, -qJD(4) * t72; t69, t69, 0, t68; t67, t67, 0, t70; 0, 0, 0, -qJD(4) * t73; -t81, -t81, 0, 0; t71, t71, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,4);
end