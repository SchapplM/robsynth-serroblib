% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% JgD_rot [3x7]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JgD_rot = S7RRRRRRR1_jacobigD_rot_5_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobigD_rot_5_floatb_twist_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_jacobigD_rot_5_floatb_twist_sym_varpar: qJD has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobigD_rot_5_floatb_twist_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-26 21:20:53
% EndTime: 2018-11-26 21:20:54
% DurationCPUTime: 0.12s
% Computational Cost: add. (34->26), mult. (110->52), div. (0->0), fcn. (112->8), ass. (0->24)
t188 = sin(qJ(4));
t193 = cos(qJ(3));
t194 = cos(qJ(2));
t199 = qJD(1) * t194 + qJD(3);
t212 = t199 * t188 * t193;
t189 = sin(qJ(3));
t200 = -qJD(3) * t194 - qJD(1);
t211 = t200 * t189;
t209 = t193 * t194;
t191 = sin(qJ(1));
t208 = qJD(1) * t191;
t195 = cos(qJ(1));
t207 = qJD(1) * t195;
t190 = sin(qJ(2));
t206 = qJD(2) * t190;
t205 = qJD(2) * t193;
t204 = qJD(2) * t194;
t203 = qJD(3) * t190;
t201 = t195 * qJD(2);
t198 = -qJD(4) + t205;
t197 = t200 * t193;
t196 = t190 * t208 - t194 * t201;
t192 = cos(qJ(4));
t1 = [0, t207, t196, t195 * t197 + (t190 * t201 + t199 * t191) * t189, -t191 * t212 + (-t198 * t190 + t211) * t188 * t195 + ((-t191 * t189 + t195 * t209) * qJD(4) + t196) * t192, 0, 0; 0, t208, -t190 * t207 - t191 * t204, t191 * t197 + (t191 * t206 - t199 * t195) * t189 (t212 + (-qJD(1) * t190 + t189 * qJD(4)) * t192) * t195 + ((-t190 * t205 + t211) * t188 - t192 * t204 + (t190 * t188 + t192 * t209) * qJD(4)) * t191, 0, 0; 0, 0, -t206, -t189 * t204 - t193 * t203 (qJD(4) * t193 - qJD(2)) * t192 * t190 + (-t189 * t203 + t198 * t194) * t188, 0, 0;];
JgD_rot  = t1;
