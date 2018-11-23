% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRR10_jacobigD_rot_3_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobigD_rot_3_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_jacobigD_rot_3_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobigD_rot_3_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:17
% EndTime: 2018-11-23 11:27:17
% DurationCPUTime: 0.04s
% Computational Cost: add. (24->14), mult. (43->29), div. (0->0), fcn. (35->11), ass. (0->15)
t191 = pkin(6) + qJ(2);
t203 = cos(t191) / 0.2e1;
t194 = sin(pkin(6));
t202 = t194 * cos(pkin(7));
t201 = qJD(1) * t194;
t200 = qJD(2) * cos(qJ(2));
t199 = cos(qJ(1));
t197 = sin(qJ(1));
t196 = sin(qJ(2));
t193 = sin(pkin(7));
t192 = pkin(6) - qJ(2);
t190 = cos(t192);
t188 = t190 / 0.2e1 + t203;
t187 = (sin(t192) / 0.2e1 - sin(t191) / 0.2e1) * qJD(2);
t1 = [0, t199 * t201 -(-t197 * t187 - t199 * t200) * t193 + (-(-t188 * t199 + t196 * t197) * t193 + t199 * t202) * qJD(1), 0, 0, 0; 0, t197 * t201 -(t199 * t187 - t197 * t200) * t193 + (-(-t188 * t197 - t196 * t199) * t193 + t197 * t202) * qJD(1), 0, 0, 0; 0, 0 -(t203 - t190 / 0.2e1) * qJD(2) * t193, 0, 0, 0;];
JgD_rot  = t1;
