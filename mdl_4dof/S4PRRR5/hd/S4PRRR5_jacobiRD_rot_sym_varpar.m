% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S4PRRR5
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% JRD_rot [9x4]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S4PRRR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_jacobiRD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_jacobiRD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_jacobiRD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:33:55
	% EndTime: 2019-12-31 16:33:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:33:55
	% EndTime: 2019-12-31 16:33:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:33:55
	% EndTime: 2019-12-31 16:33:55
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (4->4), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t20 = qJD(2) * sin(qJ(2));
	t19 = qJD(2) * cos(qJ(2));
	t16 = cos(pkin(7));
	t15 = sin(pkin(7));
	t1 = [0, -t16 * t19, 0, 0; 0, -t15 * t19, 0, 0; 0, -t20, 0, 0; 0, t16 * t20, 0, 0; 0, t15 * t20, 0, 0; 0, -t19, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:33:55
	% EndTime: 2019-12-31 16:33:55
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (32->10), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->13)
	t41 = qJ(2) + qJ(3);
	t38 = sin(t41);
	t40 = qJD(2) + qJD(3);
	t49 = t40 * t38;
	t39 = cos(t41);
	t48 = t40 * t39;
	t47 = t40 * sin(pkin(7));
	t46 = t40 * cos(pkin(7));
	t45 = t39 * t47;
	t44 = t39 * t46;
	t37 = t38 * t46;
	t36 = t38 * t47;
	t1 = [0, -t44, -t44, 0; 0, -t45, -t45, 0; 0, -t49, -t49, 0; 0, t37, t37, 0; 0, t36, t36, 0; 0, -t48, -t48, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:33:56
	% EndTime: 2019-12-31 16:33:56
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (92->19), mult. (114->36), div. (0->0), fcn. (114->6), ass. (0->29)
	t231 = qJ(2) + qJ(3);
	t228 = sin(t231);
	t230 = qJD(2) + qJD(3);
	t248 = t228 * t230;
	t234 = sin(qJ(4));
	t247 = t230 * t234;
	t235 = cos(qJ(4));
	t246 = t230 * t235;
	t232 = sin(pkin(7));
	t245 = t232 * t234;
	t244 = t232 * t235;
	t233 = cos(pkin(7));
	t243 = t233 * t234;
	t242 = t233 * t235;
	t241 = qJD(4) * t234;
	t240 = qJD(4) * t235;
	t239 = t232 * t248;
	t238 = t233 * t248;
	t229 = cos(t231);
	t237 = t228 * t240 + t229 * t247;
	t236 = t228 * t241 - t229 * t246;
	t227 = t230 * t229;
	t226 = -t228 * t246 - t229 * t241;
	t225 = t228 * t247 - t229 * t240;
	t224 = t236 * t233;
	t223 = t237 * t233;
	t222 = t236 * t232;
	t221 = t237 * t232;
	t1 = [0, t224, t224, t234 * t238 + (-t229 * t242 - t245) * qJD(4); 0, t222, t222, t234 * t239 + (-t229 * t244 + t243) * qJD(4); 0, t226, t226, -t237; 0, t223, t223, t235 * t238 + (t229 * t243 - t244) * qJD(4); 0, t221, t221, t235 * t239 + (t229 * t245 + t242) * qJD(4); 0, t225, t225, t236; 0, -t238, -t238, 0; 0, -t239, -t239, 0; 0, t227, t227, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,4);
end