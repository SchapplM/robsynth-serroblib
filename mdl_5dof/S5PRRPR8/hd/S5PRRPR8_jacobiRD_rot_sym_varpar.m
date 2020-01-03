% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRPR8
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PRRPR8_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR8_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:43:02
	% EndTime: 2019-12-31 17:43:02
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:43:02
	% EndTime: 2019-12-31 17:43:02
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:43:02
	% EndTime: 2019-12-31 17:43:02
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (4->4), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t20 = qJD(2) * sin(qJ(2));
	t19 = qJD(2) * cos(qJ(2));
	t16 = cos(pkin(8));
	t15 = sin(pkin(8));
	t1 = [0, -t16 * t19, 0, 0, 0; 0, -t15 * t19, 0, 0, 0; 0, -t20, 0, 0, 0; 0, t16 * t20, 0, 0, 0; 0, t15 * t20, 0, 0, 0; 0, -t19, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:43:02
	% EndTime: 2019-12-31 17:43:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (32->10), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->13)
	t41 = qJ(2) + qJ(3);
	t38 = sin(t41);
	t40 = qJD(2) + qJD(3);
	t49 = t40 * t38;
	t39 = cos(t41);
	t48 = t40 * t39;
	t47 = t40 * sin(pkin(8));
	t46 = t40 * cos(pkin(8));
	t45 = t39 * t47;
	t44 = t39 * t46;
	t37 = t38 * t46;
	t36 = t38 * t47;
	t1 = [0, -t44, -t44, 0, 0; 0, -t45, -t45, 0, 0; 0, -t49, -t49, 0, 0; 0, t37, t37, 0, 0; 0, t36, t36, 0, 0; 0, -t48, -t48, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:43:02
	% EndTime: 2019-12-31 17:43:02
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (44->10), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->13)
	t50 = qJ(2) + qJ(3) + pkin(9);
	t48 = sin(t50);
	t51 = qJD(2) + qJD(3);
	t59 = t51 * t48;
	t49 = cos(t50);
	t58 = t51 * t49;
	t57 = t51 * sin(pkin(8));
	t56 = t51 * cos(pkin(8));
	t55 = t49 * t57;
	t54 = t49 * t56;
	t47 = t48 * t56;
	t46 = t48 * t57;
	t1 = [0, -t54, -t54, 0, 0; 0, -t55, -t55, 0, 0; 0, -t59, -t59, 0, 0; 0, t47, t47, 0, 0; 0, t46, t46, 0, 0; 0, -t58, -t58, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:43:03
	% EndTime: 2019-12-31 17:43:03
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (134->19), mult. (114->36), div. (0->0), fcn. (114->6), ass. (0->29)
	t239 = qJ(2) + qJ(3) + pkin(9);
	t237 = sin(t239);
	t240 = qJD(2) + qJD(3);
	t257 = t237 * t240;
	t243 = sin(qJ(5));
	t256 = t240 * t243;
	t244 = cos(qJ(5));
	t255 = t240 * t244;
	t241 = sin(pkin(8));
	t254 = t241 * t243;
	t253 = t241 * t244;
	t242 = cos(pkin(8));
	t252 = t242 * t243;
	t251 = t242 * t244;
	t250 = qJD(5) * t243;
	t249 = qJD(5) * t244;
	t248 = t241 * t257;
	t247 = t242 * t257;
	t238 = cos(t239);
	t246 = t237 * t249 + t238 * t256;
	t245 = t237 * t250 - t238 * t255;
	t236 = t240 * t238;
	t235 = -t237 * t255 - t238 * t250;
	t234 = t237 * t256 - t238 * t249;
	t233 = t245 * t242;
	t232 = t246 * t242;
	t231 = t245 * t241;
	t230 = t246 * t241;
	t1 = [0, t233, t233, 0, t243 * t247 + (-t238 * t251 - t254) * qJD(5); 0, t231, t231, 0, t243 * t248 + (-t238 * t253 + t252) * qJD(5); 0, t235, t235, 0, -t246; 0, t232, t232, 0, t244 * t247 + (t238 * t252 - t253) * qJD(5); 0, t230, t230, 0, t244 * t248 + (t238 * t254 + t251) * qJD(5); 0, t234, t234, 0, t245; 0, -t247, -t247, 0, 0; 0, -t248, -t248, 0, 0; 0, t236, t236, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end