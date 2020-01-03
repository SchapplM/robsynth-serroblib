% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PPRRP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PPRRP4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRP4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_jacobiRD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:34:52
	% EndTime: 2019-12-31 17:34:52
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:34:52
	% EndTime: 2019-12-31 17:34:52
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:34:52
	% EndTime: 2019-12-31 17:34:52
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:34:52
	% EndTime: 2019-12-31 17:34:52
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
	% StartTime: 2019-12-31 17:34:52
	% EndTime: 2019-12-31 17:34:52
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (26->8), mult. (82->16), div. (0->0), fcn. (90->6), ass. (0->17)
	t68 = cos(pkin(7));
	t70 = sin(qJ(3));
	t77 = sin(pkin(7));
	t79 = cos(qJ(3));
	t64 = -t68 * t70 + t77 * t79;
	t69 = sin(qJ(4));
	t76 = qJD(4) * t69;
	t71 = cos(qJ(4));
	t75 = qJD(4) * t71;
	t61 = t64 * qJD(3);
	t63 = -t68 * t79 - t77 * t70;
	t73 = -t61 * t69 + t63 * t75;
	t72 = -t61 * t71 - t63 * t76;
	t62 = t63 * qJD(3);
	t60 = t62 * t71 - t64 * t76;
	t59 = -t62 * t69 - t64 * t75;
	t1 = [0, 0, t60, t73, 0; 0, 0, t72, t59, 0; 0, 0, 0, t76, 0; 0, 0, t59, t72, 0; 0, 0, -t73, -t60, 0; 0, 0, 0, t75, 0; 0, 0, t61, 0, 0; 0, 0, t62, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:34:53
	% EndTime: 2019-12-31 17:34:53
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (26->8), mult. (82->16), div. (0->0), fcn. (90->6), ass. (0->17)
	t78 = cos(pkin(7));
	t80 = sin(qJ(3));
	t87 = sin(pkin(7));
	t89 = cos(qJ(3));
	t74 = -t78 * t80 + t87 * t89;
	t79 = sin(qJ(4));
	t86 = qJD(4) * t79;
	t81 = cos(qJ(4));
	t85 = qJD(4) * t81;
	t71 = t74 * qJD(3);
	t73 = -t78 * t89 - t87 * t80;
	t83 = -t71 * t79 + t73 * t85;
	t82 = -t71 * t81 - t73 * t86;
	t72 = t73 * qJD(3);
	t70 = t72 * t81 - t74 * t86;
	t69 = -t72 * t79 - t74 * t85;
	t1 = [0, 0, t70, t83, 0; 0, 0, t82, t69, 0; 0, 0, 0, t86, 0; 0, 0, t69, t82, 0; 0, 0, -t83, -t70, 0; 0, 0, 0, t85, 0; 0, 0, t71, 0, 0; 0, 0, t72, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end