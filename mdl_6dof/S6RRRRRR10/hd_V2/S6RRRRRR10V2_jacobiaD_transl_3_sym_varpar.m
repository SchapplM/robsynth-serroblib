% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR10V2_jacobiaD_transl_3_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_transl_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_transl_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR10V2_jacobiaD_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_transl_3_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:56:31
% EndTime: 2019-04-11 14:56:32
% DurationCPUTime: 0.09s
% Computational Cost: add. (77->25), mult. (110->37), div. (0->0), fcn. (71->6), ass. (0->26)
t36 = qJD(2) + qJD(3);
t37 = qJ(2) + qJ(3);
t35 = cos(t37);
t55 = r_i_i_C(2) * t35;
t34 = sin(t37);
t57 = r_i_i_C(1) * t34;
t46 = t55 + t57;
t44 = t46 * t36;
t38 = sin(qJ(2));
t58 = pkin(2) * t38;
t59 = qJD(2) * t58 + t44;
t56 = r_i_i_C(2) * t34;
t54 = t35 * t36;
t39 = sin(qJ(1));
t53 = qJD(1) * t39;
t41 = cos(qJ(1));
t52 = qJD(1) * t41;
t40 = cos(qJ(2));
t51 = qJD(2) * t40;
t50 = r_i_i_C(1) * t54;
t49 = t36 * t56;
t48 = qJD(1) * t55;
t45 = -t40 * pkin(2) - r_i_i_C(1) * t35 - pkin(1) + t56;
t43 = t39 * t48 + t53 * t57 + (t49 - t50) * t41;
t29 = t39 * t49;
t1 = [t59 * t39 + (-r_i_i_C(3) * t39 + t45 * t41) * qJD(1) (t38 * t53 - t41 * t51) * pkin(2) + t43, t43, 0, 0, 0; -t59 * t41 + (r_i_i_C(3) * t41 + t45 * t39) * qJD(1), t29 + (-pkin(2) * t51 - t50) * t39 + (-t46 - t58) * t52, -t41 * t48 + t29 + (-t34 * t52 - t39 * t54) * r_i_i_C(1), 0, 0, 0; 0, -t59, -t44, 0, 0, 0;];
JaD_transl  = t1;
