% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRPRR9
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRPRR9_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRPRR9_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR9_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:39:59
	% EndTime: 2019-12-31 17:39:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:39:59
	% EndTime: 2019-12-31 17:39:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:39:59
	% EndTime: 2019-12-31 17:39:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->3), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->4)
	t32 = pkin(8) + qJ(2);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [0, (-r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(2), 0, 0, 0; 0, (-r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(2), 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:39:59
	% EndTime: 2019-12-31 17:39:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (20->7), mult. (20->10), div. (0->0), fcn. (12->2), ass. (0->6)
	t13 = -pkin(2) - r_i_i_C(1);
	t12 = r_i_i_C(3) + qJ(3);
	t11 = pkin(8) + qJ(2);
	t10 = cos(t11);
	t9 = sin(t11);
	t1 = [0, t10 * qJD(3) + (t13 * t10 - t12 * t9) * qJD(2), qJD(2) * t10, 0, 0; 0, t9 * qJD(3) + (t12 * t10 + t13 * t9) * qJD(2), qJD(2) * t9, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:39:59
	% EndTime: 2019-12-31 17:39:59
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (78->16), mult. (88->25), div. (0->0), fcn. (74->4), ass. (0->15)
	t92 = -pkin(2) - pkin(3);
	t83 = pkin(8) + qJ(2);
	t81 = sin(t83);
	t91 = qJD(2) * t81;
	t82 = cos(t83);
	t90 = qJD(2) * t82;
	t84 = sin(qJ(4));
	t89 = qJD(4) * t84;
	t85 = cos(qJ(4));
	t88 = qJD(4) * t85;
	t75 = -t81 * t89 - t82 * t88 + (t81 * t84 + t82 * t85) * qJD(2);
	t76 = t81 * t88 - t82 * t89 + t84 * t90 - t85 * t91;
	t87 = -t76 * r_i_i_C(1) - t75 * r_i_i_C(2);
	t86 = t75 * r_i_i_C(1) - t76 * r_i_i_C(2);
	t1 = [0, t82 * qJD(3) + (-qJ(3) * t81 + t92 * t82) * qJD(2) - t86, t90, t86, 0; 0, t81 * qJD(3) + (qJ(3) * t82 + t92 * t81) * qJD(2) - t87, t91, t87, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:39:59
	% EndTime: 2019-12-31 17:39:59
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (253->24), mult. (308->39), div. (0->0), fcn. (276->6), ass. (0->22)
	t85 = qJD(2) - qJD(4);
	t72 = pkin(8) + qJ(2);
	t70 = sin(t72);
	t71 = cos(t72);
	t75 = sin(qJ(4));
	t76 = cos(qJ(4));
	t51 = -t70 * t75 - t71 * t76;
	t49 = t85 * t51;
	t52 = -t70 * t76 + t71 * t75;
	t50 = t85 * t52;
	t60 = sin(qJ(5));
	t61 = cos(qJ(5));
	t68 = r_i_i_C(1) * t61 - r_i_i_C(2) * t60 + pkin(4);
	t77 = pkin(7) + r_i_i_C(3);
	t82 = (r_i_i_C(1) * t60 + r_i_i_C(2) * t61) * qJD(5);
	t84 = -t77 * t49 - t68 * t50 - t51 * t82;
	t83 = -t68 * t49 + t77 * t50 + t52 * t82;
	t81 = qJ(3) * qJD(2);
	t78 = qJD(3) + (-pkin(2) - pkin(3)) * qJD(2);
	t74 = qJD(5) * t60;
	t73 = qJD(5) * t61;
	t1 = [0, -t70 * t81 + t78 * t71 - t83, qJD(2) * t71, t83, (-t50 * t61 - t51 * t74) * r_i_i_C(2) + (-t50 * t60 + t51 * t73) * r_i_i_C(1); 0, t78 * t70 + t71 * t81 - t84, qJD(2) * t70, t84, (t49 * t61 - t52 * t74) * r_i_i_C(2) + (t49 * t60 + t52 * t73) * r_i_i_C(1); 0, 0, 0, 0, t82;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end