% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRPR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:42
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPRPR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:35
	% EndTime: 2019-10-24 10:42:35
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:35
	% EndTime: 2019-10-24 10:42:35
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = qJD(1) * cos(qJ(1));
	t7 = qJD(1) * sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; t7, 0, 0, 0, 0; -t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; t9, 0, 0, 0, 0; t7, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:35
	% EndTime: 2019-10-24 10:42:35
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (5->2), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = qJ(1) + pkin(8);
	t13 = qJD(1) * cos(t12);
	t10 = qJD(1) * sin(t12);
	t1 = [0, 0, 0, 0, 0; t10, 0, 0, 0, 0; -t13, 0, 0, 0, 0; 0, 0, 0, 0, 0; t13, 0, 0, 0, 0; t10, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:35
	% EndTime: 2019-10-24 10:42:35
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (30->11), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t85 = sin(qJ(3));
	t90 = qJD(1) * t85;
	t86 = cos(qJ(3));
	t89 = qJD(1) * t86;
	t88 = qJD(3) * t85;
	t87 = qJD(3) * t86;
	t84 = qJ(1) + pkin(8);
	t83 = cos(t84);
	t82 = sin(t84);
	t81 = -t82 * t88 + t83 * t89;
	t80 = t82 * t87 + t83 * t90;
	t79 = t82 * t89 + t83 * t88;
	t78 = t82 * t90 - t83 * t87;
	t1 = [0, 0, -t88, 0, 0; t79, 0, t80, 0, 0; -t81, 0, t78, 0, 0; 0, 0, -t87, 0, 0; -t78, 0, t81, 0, 0; t80, 0, t79, 0, 0; 0, 0, 0, 0, 0; -qJD(1) * t83, 0, 0, 0, 0; -qJD(1) * t82, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:35
	% EndTime: 2019-10-24 10:42:35
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (48->12), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->16)
	t102 = qJ(1) + pkin(8);
	t98 = sin(t102);
	t107 = qJD(1) * t98;
	t101 = qJ(3) + pkin(9);
	t97 = sin(t101);
	t106 = qJD(3) * t97;
	t99 = cos(t101);
	t105 = qJD(3) * t99;
	t100 = cos(t102);
	t104 = qJD(1) * t100;
	t103 = qJD(3) * t100;
	t96 = t99 * t104 - t98 * t106;
	t95 = t97 * t104 + t98 * t105;
	t94 = t97 * t103 + t99 * t107;
	t93 = -t99 * t103 + t97 * t107;
	t1 = [0, 0, -t106, 0, 0; t94, 0, t95, 0, 0; -t96, 0, t93, 0, 0; 0, 0, -t105, 0, 0; -t93, 0, t96, 0, 0; t95, 0, t94, 0, 0; 0, 0, 0, 0, 0; -t104, 0, 0, 0, 0; -t107, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:35
	% EndTime: 2019-10-24 10:42:35
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (116->15), mult. (54->12), div. (0->0), fcn. (54->4), ass. (0->16)
	t145 = qJ(3) + pkin(9) + qJ(5);
	t141 = sin(t145);
	t146 = qJD(3) + qJD(5);
	t151 = t146 * t141;
	t142 = cos(t145);
	t150 = t146 * t142;
	t147 = qJ(1) + pkin(8);
	t143 = sin(t147);
	t149 = qJD(1) * t143;
	t144 = cos(t147);
	t148 = qJD(1) * t144;
	t140 = t142 * t148 - t143 * t151;
	t139 = t141 * t148 + t143 * t150;
	t138 = t142 * t149 + t144 * t151;
	t137 = t141 * t149 - t144 * t150;
	t1 = [0, 0, -t151, 0, -t151; t138, 0, t139, 0, t139; -t140, 0, t137, 0, t137; 0, 0, -t150, 0, -t150; -t137, 0, t140, 0, t140; t139, 0, t138, 0, t138; 0, 0, 0, 0, 0; -t148, 0, 0, 0, 0; -t149, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end