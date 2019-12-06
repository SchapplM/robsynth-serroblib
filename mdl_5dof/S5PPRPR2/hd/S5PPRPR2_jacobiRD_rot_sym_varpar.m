% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PPRPR2
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PPRPR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRPR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:03:38
	% EndTime: 2019-12-05 15:03:38
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:03:38
	% EndTime: 2019-12-05 15:03:38
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:03:38
	% EndTime: 2019-12-05 15:03:38
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:03:38
	% EndTime: 2019-12-05 15:03:38
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
	% StartTime: 2019-12-05 15:03:39
	% EndTime: 2019-12-05 15:03:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->3), mult. (10->8), div. (0->0), fcn. (10->4), ass. (0->6)
	t109 = qJD(3) * sin(pkin(7));
	t108 = qJD(3) * cos(pkin(7));
	t105 = pkin(8) + qJ(3);
	t104 = cos(t105);
	t103 = sin(t105);
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, t104 * t108, 0, 0; 0, 0, t104 * t109, 0, 0; 0, 0, qJD(3) * t103, 0, 0; 0, 0, -t103 * t108, 0, 0; 0, 0, -t103 * t109, 0, 0; 0, 0, qJD(3) * t104, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:03:39
	% EndTime: 2019-12-05 15:03:39
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (46->18), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->21)
	t158 = sin(pkin(7));
	t160 = sin(qJ(5));
	t174 = t158 * t160;
	t161 = cos(qJ(5));
	t173 = t158 * t161;
	t159 = cos(pkin(7));
	t172 = t159 * t160;
	t171 = t159 * t161;
	t157 = pkin(8) + qJ(3);
	t156 = cos(t157);
	t170 = qJD(3) * t156;
	t169 = qJD(3) * t160;
	t168 = qJD(3) * t161;
	t167 = qJD(5) * t160;
	t166 = qJD(5) * t161;
	t165 = t158 * t170;
	t164 = t159 * t170;
	t155 = sin(t157);
	t163 = t155 * t168 + t156 * t167;
	t162 = -t155 * t169 + t156 * t166;
	t1 = [0, 0, t162 * t159, 0, t161 * t164 + (-t155 * t172 - t173) * qJD(5); 0, 0, t162 * t158, 0, t161 * t165 + (-t155 * t174 + t171) * qJD(5); 0, 0, t155 * t166 + t156 * t169, 0, t163; 0, 0, -t163 * t159, 0, -t160 * t164 + (-t155 * t171 + t174) * qJD(5); 0, 0, -t163 * t158, 0, -t160 * t165 + (-t155 * t173 - t172) * qJD(5); 0, 0, -t155 * t167 + t156 * t168, 0, t162; 0, 0, -t164, 0, 0; 0, 0, -t165, 0, 0; 0, 0, -qJD(3) * t155, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end