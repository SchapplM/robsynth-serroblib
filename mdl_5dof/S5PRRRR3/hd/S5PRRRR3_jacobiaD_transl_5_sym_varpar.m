% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
%
% Output:
% JaD_transl [3x5]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-06-06 15:46
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRRR3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_jacobiaD_transl_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_jacobiaD_transl_5_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRR3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR3_jacobiaD_transl_5_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-06-06 15:46:08
% EndTime: 2019-06-06 15:46:08
% DurationCPUTime: 0.08s
% Computational Cost: add. (195->21), mult. (136->35), div. (0->0), fcn. (84->8), ass. (0->24)
t71 = pkin(6) + r_i_i_C(3);
t49 = qJ(2) + qJ(3);
t47 = qJ(4) + t49;
t42 = sin(t47);
t43 = cos(t47);
t51 = cos(qJ(5));
t62 = qJD(5) * t51;
t48 = qJD(2) + qJD(3);
t44 = qJD(4) + t48;
t50 = sin(qJ(5));
t66 = t44 * t50;
t70 = t42 * t66 - t43 * t62;
t69 = t42 * t62 + t43 * t66;
t68 = pkin(3) * t48;
t65 = t44 * t51;
t64 = pkin(2) * qJD(2);
t63 = qJD(5) * t50;
t59 = t42 * t63;
t56 = t42 * t65 + t43 * t63;
t55 = r_i_i_C(1) * t59 + (-r_i_i_C(1) * t43 * t51 - t71 * t42) * t44 + t69 * r_i_i_C(2);
t54 = -cos(t49) * t68 + t55;
t53 = t71 * t43 * t44 - t56 * r_i_i_C(1) + t70 * r_i_i_C(2);
t52 = -sin(t49) * t68 + t53;
t1 = [0, -cos(qJ(2)) * t64 + t54, t54, t55, t70 * r_i_i_C(1) + t56 * r_i_i_C(2); 0, -sin(qJ(2)) * t64 + t52, t52, t53, (-t43 * t65 + t59) * r_i_i_C(2) - t69 * r_i_i_C(1); 0, 0, 0, 0, (-r_i_i_C(1) * t50 - r_i_i_C(2) * t51) * qJD(5);];
JaD_transl  = t1;
