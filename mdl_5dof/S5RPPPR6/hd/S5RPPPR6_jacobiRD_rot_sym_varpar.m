% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPPPR6_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR6_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:48:16
	% EndTime: 2019-12-31 17:48:16
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:48:16
	% EndTime: 2019-12-31 17:48:16
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:48:16
	% EndTime: 2019-12-31 17:48:16
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (3->3), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t22 = qJD(1) * sin(qJ(1));
	t21 = qJD(1) * cos(qJ(1));
	t18 = cos(pkin(7));
	t17 = sin(pkin(7));
	t1 = [-t18 * t21, 0, 0, 0, 0; -t18 * t22, 0, 0, 0, 0; 0, 0, 0, 0, 0; t17 * t21, 0, 0, 0, 0; t17 * t22, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t22, 0, 0, 0, 0; t21, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:48:16
	% EndTime: 2019-12-31 17:48:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t104 = qJD(1) * sin(qJ(1));
	t103 = qJD(1) * cos(qJ(1));
	t100 = cos(pkin(7));
	t99 = sin(pkin(7));
	t1 = [-t104, 0, 0, 0, 0; t103, 0, 0, 0, 0; 0, 0, 0, 0, 0; t100 * t103, 0, 0, 0, 0; t100 * t104, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t99 * t103, 0, 0, 0, 0; -t99 * t104, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:48:16
	% EndTime: 2019-12-31 17:48:16
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->6), mult. (24->15), div. (0->0), fcn. (24->6), ass. (0->11)
	t124 = sin(pkin(8));
	t128 = sin(qJ(1));
	t134 = t124 * t128;
	t129 = cos(qJ(1));
	t133 = t124 * t129;
	t126 = cos(pkin(8));
	t132 = t126 * t128;
	t131 = t126 * t129;
	t130 = qJD(1) * cos(pkin(7));
	t125 = sin(pkin(7));
	t1 = [(-t125 * t133 - t132) * qJD(1), 0, 0, 0, 0; (-t125 * t134 + t131) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0; (-t125 * t131 + t134) * qJD(1), 0, 0, 0, 0; (-t125 * t132 - t133) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0; -t129 * t130, 0, 0, 0, 0; -t128 * t130, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:48:17
	% EndTime: 2019-12-31 17:48:17
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->20), mult. (156->44), div. (0->0), fcn. (164->8), ass. (0->27)
	t236 = sin(pkin(7));
	t237 = cos(pkin(8));
	t240 = sin(qJ(1));
	t250 = t240 * t237;
	t235 = sin(pkin(8));
	t242 = cos(qJ(1));
	t254 = t235 * t242;
	t232 = t236 * t254 + t250;
	t230 = t232 * qJD(1);
	t249 = t242 * t237;
	t255 = t235 * t240;
	t233 = -t236 * t255 + t249;
	t239 = sin(qJ(5));
	t241 = cos(qJ(5));
	t238 = cos(pkin(7));
	t248 = qJD(1) * t238;
	t246 = t242 * t248;
	t253 = t238 * t239;
	t256 = (-t233 * t241 + t240 * t253) * qJD(5) + t230 * t239 - t241 * t246;
	t252 = t238 * t241;
	t251 = t238 * t242;
	t247 = t240 * t248;
	t243 = -t239 * t246 - t230 * t241 + (-t233 * t239 - t240 * t252) * qJD(5);
	t231 = t233 * qJD(1);
	t229 = -t239 * t247 + t231 * t241 + (-t232 * t239 + t241 * t251) * qJD(5);
	t228 = -t241 * t247 - t231 * t239 + (-t232 * t241 - t239 * t251) * qJD(5);
	t1 = [t243, 0, 0, 0, t228; t229, 0, 0, 0, -t256; 0, 0, 0, 0, (t235 * t252 - t236 * t239) * qJD(5); t256, 0, 0, 0, -t229; t228, 0, 0, 0, t243; 0, 0, 0, 0, (-t235 * t253 - t236 * t241) * qJD(5); (t236 * t249 - t255) * qJD(1), 0, 0, 0, 0; (t236 * t250 + t254) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end