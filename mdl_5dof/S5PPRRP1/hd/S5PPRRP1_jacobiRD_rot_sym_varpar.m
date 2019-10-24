% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PPRRP1
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:19
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PPRRP1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRP1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:18
	% EndTime: 2019-10-24 10:19:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:18
	% EndTime: 2019-10-24 10:19:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:18
	% EndTime: 2019-10-24 10:19:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:18
	% EndTime: 2019-10-24 10:19:18
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
	% StartTime: 2019-10-24 10:19:19
	% EndTime: 2019-10-24 10:19:19
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (45->16), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->21)
	t166 = sin(pkin(7));
	t168 = sin(qJ(4));
	t182 = t166 * t168;
	t169 = cos(qJ(4));
	t181 = t166 * t169;
	t167 = cos(pkin(7));
	t180 = t167 * t168;
	t179 = t167 * t169;
	t165 = pkin(8) + qJ(3);
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
	t1 = [0, 0, t171 * t167, t168 * t172 + (-t164 * t179 - t182) * qJD(4), 0; 0, 0, t171 * t166, t168 * t173 + (-t164 * t181 + t180) * qJD(4), 0; 0, 0, -t163 * t176 - t164 * t175, -t170, 0; 0, 0, t170 * t167, t169 * t172 + (t164 * t180 - t181) * qJD(4), 0; 0, 0, t170 * t166, t169 * t173 + (t164 * t182 + t179) * qJD(4), 0; 0, 0, t163 * t177 - t164 * t174, t171, 0; 0, 0, -t172, 0, 0; 0, 0, -t173, 0, 0; 0, 0, qJD(3) * t164, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:19
	% EndTime: 2019-10-24 10:19:19
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (45->16), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->21)
	t179 = sin(pkin(7));
	t181 = sin(qJ(4));
	t195 = t179 * t181;
	t182 = cos(qJ(4));
	t194 = t179 * t182;
	t180 = cos(pkin(7));
	t193 = t180 * t181;
	t192 = t180 * t182;
	t178 = pkin(8) + qJ(3);
	t176 = sin(t178);
	t191 = qJD(3) * t176;
	t190 = qJD(3) * t181;
	t189 = qJD(3) * t182;
	t188 = qJD(4) * t181;
	t187 = qJD(4) * t182;
	t186 = t179 * t191;
	t185 = t180 * t191;
	t177 = cos(t178);
	t184 = t176 * t188 - t177 * t189;
	t183 = t176 * t187 + t177 * t190;
	t1 = [0, 0, t184 * t180, t181 * t185 + (-t177 * t192 - t195) * qJD(4), 0; 0, 0, t184 * t179, t181 * t186 + (-t177 * t194 + t193) * qJD(4), 0; 0, 0, -t176 * t189 - t177 * t188, -t183, 0; 0, 0, t183 * t180, t182 * t185 + (t177 * t193 - t194) * qJD(4), 0; 0, 0, t183 * t179, t182 * t186 + (t177 * t195 + t192) * qJD(4), 0; 0, 0, t176 * t190 - t177 * t187, t184, 0; 0, 0, -t185, 0, 0; 0, 0, -t186, 0, 0; 0, 0, qJD(3) * t177, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end