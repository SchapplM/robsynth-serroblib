% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4PPRR1
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% JaD_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:23
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4PPRR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_jacobiaD_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_jacobiaD_transl_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4PPRR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PPRR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_jacobiaD_transl_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:23:31
	% EndTime: 2019-10-09 20:23:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:23:31
	% EndTime: 2019-10-09 20:23:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:23:31
	% EndTime: 2019-10-09 20:23:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:23:31
	% EndTime: 2019-10-09 20:23:31
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (6->4), mult. (20->10), div. (0->0), fcn. (16->4), ass. (0->7)
	t49 = cos(qJ(3));
	t48 = sin(qJ(3));
	t47 = cos(pkin(6));
	t46 = sin(pkin(6));
	t45 = (t46 * t48 + t47 * t49) * qJD(3);
	t44 = (-t46 * t49 + t47 * t48) * qJD(3);
	t1 = [0, 0, -t45 * r_i_i_C(1) + t44 * r_i_i_C(2), 0; 0, 0, t44 * r_i_i_C(1) + t45 * r_i_i_C(2), 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:23:31
	% EndTime: 2019-10-09 20:23:31
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (48->10), mult. (52->17), div. (0->0), fcn. (40->6), ass. (0->14)
	t68 = qJ(3) + qJ(4);
	t65 = sin(t68);
	t66 = cos(t68);
	t67 = qJD(3) + qJD(4);
	t69 = sin(pkin(6));
	t70 = cos(pkin(6));
	t63 = (t65 * t70 - t66 * t69) * t67;
	t64 = (t65 * t69 + t66 * t70) * t67;
	t75 = t63 * r_i_i_C(1) + t64 * r_i_i_C(2);
	t74 = -t64 * r_i_i_C(1) + t63 * r_i_i_C(2);
	t73 = pkin(3) * qJD(3);
	t72 = cos(qJ(3));
	t71 = sin(qJ(3));
	t1 = [0, 0, (-t69 * t71 - t70 * t72) * t73 + t74, t74; 0, 0, (-t69 * t72 + t70 * t71) * t73 + t75, t75; 0, 0, 0, 0;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,4);
end