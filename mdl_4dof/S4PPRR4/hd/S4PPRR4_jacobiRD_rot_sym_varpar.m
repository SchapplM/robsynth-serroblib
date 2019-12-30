% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S4PPRR4
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
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
% 
% Output:
% JRD_rot [9x4]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 11:56
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S4PPRR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_jacobiRD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_jacobiRD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PPRR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_jacobiRD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 11:56:34
	% EndTime: 2019-12-29 11:56:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 11:56:34
	% EndTime: 2019-12-29 11:56:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 11:56:29
	% EndTime: 2019-12-29 11:56:29
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 11:56:34
	% EndTime: 2019-12-29 11:56:34
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (10->5), mult. (10->8), div. (0->0), fcn. (10->4), ass. (0->6)
	t24 = qJD(3) * sin(pkin(6));
	t23 = qJD(3) * cos(pkin(6));
	t20 = pkin(7) + qJ(3);
	t19 = cos(t20);
	t18 = sin(t20);
	t1 = [0, 0, -t19 * t23, 0; 0, 0, -t19 * t24, 0; 0, 0, -qJD(3) * t18, 0; 0, 0, t18 * t23, 0; 0, 0, t18 * t24, 0; 0, 0, -qJD(3) * t19, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 11:56:35
	% EndTime: 2019-12-29 11:56:36
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (45->16), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->21)
	t166 = sin(pkin(6));
	t168 = sin(qJ(4));
	t182 = t166 * t168;
	t169 = cos(qJ(4));
	t181 = t166 * t169;
	t167 = cos(pkin(6));
	t180 = t167 * t168;
	t179 = t167 * t169;
	t165 = pkin(7) + qJ(3);
	t163 = sin(t165);
	t178 = qJD(3) * t163;
	t177 = qJD(3) * t168;
	t176 = qJD(3) * t169;
	t175 = qJD(4) * t168;
	t174 = qJD(4) * t169;
	t173 = t166 * t178;
	t172 = t167 * t178;
	t164 = cos(t165);
	t171 = t163 * t175 - t164 * t176;
	t170 = t163 * t174 + t164 * t177;
	t1 = [0, 0, t171 * t167, t168 * t172 + (-t164 * t179 - t182) * qJD(4); 0, 0, t171 * t166, t168 * t173 + (-t164 * t181 + t180) * qJD(4); 0, 0, -t163 * t176 - t164 * t175, -t170; 0, 0, t170 * t167, t169 * t172 + (t164 * t180 - t181) * qJD(4); 0, 0, t170 * t166, t169 * t173 + (t164 * t182 + t179) * qJD(4); 0, 0, t163 * t177 - t164 * t174, t171; 0, 0, -t172, 0; 0, 0, -t173, 0; 0, 0, qJD(3) * t164, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,4);
end