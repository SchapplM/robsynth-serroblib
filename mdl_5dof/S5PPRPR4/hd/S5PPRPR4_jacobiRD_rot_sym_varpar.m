% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PPRPR4
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PPRPR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRPR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:32:39
	% EndTime: 2019-12-31 17:32:39
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:32:39
	% EndTime: 2019-12-31 17:32:39
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:32:39
	% EndTime: 2019-12-31 17:32:39
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:32:39
	% EndTime: 2019-12-31 17:32:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->3), mult. (16->6), div. (0->0), fcn. (16->4), ass. (0->7)
	t55 = cos(qJ(3));
	t54 = sin(qJ(3));
	t53 = cos(pkin(7));
	t52 = sin(pkin(7));
	t51 = (t52 * t54 + t53 * t55) * qJD(3);
	t50 = (-t52 * t55 + t53 * t54) * qJD(3);
	t1 = [0, 0, -t51, 0, 0; 0, 0, t50, 0, 0; 0, 0, 0, 0, 0; 0, 0, t50, 0, 0; 0, 0, t51, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:32:39
	% EndTime: 2019-12-31 17:32:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (28->10), div. (0->0), fcn. (28->6), ass. (0->9)
	t44 = cos(qJ(3));
	t43 = sin(qJ(3));
	t42 = cos(pkin(7));
	t41 = cos(pkin(8));
	t40 = sin(pkin(7));
	t39 = sin(pkin(8));
	t38 = (-t40 * t43 - t42 * t44) * qJD(3);
	t37 = (-t40 * t44 + t42 * t43) * qJD(3);
	t1 = [0, 0, t38 * t41, 0, 0; 0, 0, t37 * t41, 0, 0; 0, 0, 0, 0, 0; 0, 0, -t38 * t39, 0, 0; 0, 0, -t37 * t39, 0, 0; 0, 0, 0, 0, 0; 0, 0, -t37, 0, 0; 0, 0, t38, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:32:39
	% EndTime: 2019-12-31 17:32:39
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (44->9), mult. (82->16), div. (0->0), fcn. (90->6), ass. (0->18)
	t79 = cos(pkin(7));
	t80 = sin(qJ(3));
	t86 = sin(pkin(7));
	t88 = cos(qJ(3));
	t72 = -t79 * t80 + t86 * t88;
	t78 = pkin(8) + qJ(5);
	t76 = sin(t78);
	t85 = qJD(5) * t76;
	t77 = cos(t78);
	t84 = qJD(5) * t77;
	t69 = t72 * qJD(3);
	t71 = -t79 * t88 - t86 * t80;
	t82 = -t69 * t76 + t71 * t84;
	t81 = -t69 * t77 - t71 * t85;
	t70 = t71 * qJD(3);
	t68 = t70 * t77 - t72 * t85;
	t67 = -t70 * t76 - t72 * t84;
	t1 = [0, 0, t68, 0, t82; 0, 0, t81, 0, t67; 0, 0, 0, 0, t85; 0, 0, t67, 0, t81; 0, 0, -t82, 0, -t68; 0, 0, 0, 0, t84; 0, 0, t69, 0, 0; 0, 0, t70, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end