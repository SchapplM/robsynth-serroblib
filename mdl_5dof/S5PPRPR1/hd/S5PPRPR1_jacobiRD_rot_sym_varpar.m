% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PPRPR1
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PPRPR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRPR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:01:48
	% EndTime: 2019-12-05 15:01:48
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:01:48
	% EndTime: 2019-12-05 15:01:48
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:01:48
	% EndTime: 2019-12-05 15:01:48
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:01:48
	% EndTime: 2019-12-05 15:01:48
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (10->5), mult. (10->8), div. (0->0), fcn. (10->4), ass. (0->6)
	t24 = qJD(3) * sin(pkin(7));
	t23 = qJD(3) * cos(pkin(7));
	t20 = pkin(8) + qJ(3);
	t19 = cos(t20);
	t18 = sin(t20);
	t1 = [0, 0, -t19 * t23, 0, 0; 0, 0, -t19 * t24, 0, 0; 0, 0, -qJD(3) * t18, 0, 0; 0, 0, t18 * t23, 0, 0; 0, 0, t18 * t24, 0, 0; 0, 0, -qJD(3) * t19, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:01:49
	% EndTime: 2019-12-05 15:01:49
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (14->6), mult. (21->14), div. (0->0), fcn. (21->6), ass. (0->11)
	t131 = pkin(8) + qJ(3);
	t129 = sin(t131);
	t140 = qJD(3) * t129;
	t139 = qJD(3) * sin(pkin(7));
	t138 = qJD(3) * cos(pkin(7));
	t130 = cos(t131);
	t137 = t130 * t139;
	t136 = t130 * t138;
	t134 = cos(pkin(9));
	t132 = sin(pkin(9));
	t1 = [0, 0, -t134 * t136, 0, 0; 0, 0, -t134 * t137, 0, 0; 0, 0, -t134 * t140, 0, 0; 0, 0, t132 * t136, 0, 0; 0, 0, t132 * t137, 0, 0; 0, 0, t132 * t140, 0, 0; 0, 0, -t129 * t138, 0, 0; 0, 0, -t129 * t139, 0, 0; 0, 0, qJD(3) * t130, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:01:49
	% EndTime: 2019-12-05 15:01:49
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (73->17), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->19)
	t177 = pkin(8) + qJ(3);
	t175 = cos(t177);
	t178 = sin(pkin(7));
	t189 = t175 * t178;
	t179 = cos(pkin(7));
	t188 = t175 * t179;
	t173 = sin(t177);
	t187 = qJD(3) * t173;
	t186 = qJD(3) * t175;
	t185 = qJD(5) * t173;
	t184 = qJD(5) * t175;
	t183 = t178 * t187;
	t182 = t179 * t187;
	t176 = pkin(9) + qJ(5);
	t172 = sin(t176);
	t174 = cos(t176);
	t181 = t172 * t185 - t174 * t186;
	t180 = t172 * t186 + t174 * t185;
	t1 = [0, 0, t181 * t179, 0, t172 * t182 + (-t172 * t178 - t174 * t188) * qJD(5); 0, 0, t181 * t178, 0, t172 * t183 + (t172 * t179 - t174 * t189) * qJD(5); 0, 0, -t172 * t184 - t174 * t187, 0, -t180; 0, 0, t180 * t179, 0, t174 * t182 + (t172 * t188 - t174 * t178) * qJD(5); 0, 0, t180 * t178, 0, t174 * t183 + (t172 * t189 + t174 * t179) * qJD(5); 0, 0, t172 * t187 - t174 * t184, 0, t181; 0, 0, -t182, 0, 0; 0, 0, -t183, 0, 0; 0, 0, t186, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end