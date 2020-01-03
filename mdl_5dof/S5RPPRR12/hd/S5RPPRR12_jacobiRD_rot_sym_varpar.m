% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPRR12
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPPRR12_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR12_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:07:49
	% EndTime: 2019-12-31 18:07:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:07:49
	% EndTime: 2019-12-31 18:07:49
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
	% StartTime: 2019-12-31 18:07:49
	% EndTime: 2019-12-31 18:07:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t13 = qJD(1) * sin(qJ(1));
	t11 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; t11, 0, 0, 0, 0; t13, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t13, 0, 0, 0, 0; t11, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:07:49
	% EndTime: 2019-12-31 18:07:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (4->4), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t24 = qJD(1) * sin(qJ(1));
	t23 = qJD(1) * cos(qJ(1));
	t20 = cos(pkin(8));
	t19 = sin(pkin(8));
	t1 = [-t19 * t24, 0, 0, 0, 0; t19 * t23, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t20 * t24, 0, 0, 0, 0; t20 * t23, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t23, 0, 0, 0, 0; -t24, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:07:49
	% EndTime: 2019-12-31 18:07:49
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t44 = sin(qJ(1));
	t49 = qJD(1) * t44;
	t45 = cos(qJ(1));
	t48 = qJD(1) * t45;
	t47 = qJD(4) * t44;
	t46 = qJD(4) * t45;
	t43 = pkin(8) + qJ(4);
	t42 = cos(t43);
	t41 = sin(t43);
	t40 = -t41 * t47 + t42 * t48;
	t39 = t41 * t48 + t42 * t47;
	t38 = t41 * t46 + t42 * t49;
	t37 = -t41 * t49 + t42 * t46;
	t1 = [t37, 0, 0, t40, 0; t39, 0, 0, t38, 0; 0, 0, 0, -qJD(4) * t42, 0; -t38, 0, 0, -t39, 0; t40, 0, 0, t37, 0; 0, 0, 0, qJD(4) * t41, 0; -t48, 0, 0, 0, 0; -t49, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:07:50
	% EndTime: 2019-12-31 18:07:50
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (102->29), mult. (173->57), div. (0->0), fcn. (173->6), ass. (0->31)
	t249 = cos(qJ(1));
	t245 = pkin(8) + qJ(4);
	t243 = sin(t245);
	t251 = qJD(1) * t243 + qJD(5);
	t268 = t251 * t249;
	t247 = sin(qJ(1));
	t244 = cos(t245);
	t261 = qJD(4) * t249;
	t253 = t244 * t261;
	t267 = t251 * t247 - t253;
	t266 = qJD(1) * t247;
	t265 = qJD(1) * t249;
	t246 = sin(qJ(5));
	t264 = qJD(4) * t246;
	t263 = qJD(4) * t247;
	t248 = cos(qJ(5));
	t262 = qJD(4) * t248;
	t260 = qJD(5) * t246;
	t259 = qJD(5) * t248;
	t258 = qJD(5) * t249;
	t257 = t244 * t262;
	t256 = t243 * t263;
	t255 = t244 * t263;
	t254 = t243 * t261;
	t252 = -qJD(5) * t243 - qJD(1);
	t250 = t252 * t249;
	t242 = t248 * t268 + (t252 * t246 + t257) * t247;
	t241 = t252 * t248 * t247 + (-t255 - t268) * t246;
	t240 = t246 * t250 - t267 * t248;
	t239 = t267 * t246 + t248 * t250;
	t1 = [t240, 0, 0, -t248 * t256 + (-t247 * t260 + t248 * t265) * t244, t241; t242, 0, 0, t248 * t254 + (t246 * t258 + t248 * t266) * t244, -t239; 0, 0, 0, t243 * t260 - t257, t243 * t264 - t244 * t259; t239, 0, 0, t246 * t256 + (-t246 * t265 - t247 * t259) * t244, -t242; t241, 0, 0, -t246 * t254 + (-t246 * t266 + t248 * t258) * t244, t240; 0, 0, 0, t243 * t259 + t244 * t264, t243 * t262 + t244 * t260; t244 * t266 + t254, 0, 0, t243 * t265 + t255, 0; -t244 * t265 + t256, 0, 0, t243 * t266 - t253, 0; 0, 0, 0, -qJD(4) * t243, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end