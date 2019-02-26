% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRRR10_jacobig_rot_3_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobig_rot_3_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobig_rot_3_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:17
% EndTime: 2018-11-23 11:27:17
% DurationCPUTime: 0.03s
% Computational Cost: add. (15->10), mult. (20->16), div. (0->0), fcn. (28->11), ass. (0->13)
t129 = sin(pkin(6));
t133 = sin(qJ(1));
t136 = t133 * t129;
t134 = cos(qJ(1));
t135 = t134 * t129;
t132 = sin(qJ(2));
t131 = cos(pkin(6));
t130 = cos(pkin(7));
t128 = sin(pkin(7));
t127 = pkin(6) - qJ(2);
t126 = pkin(6) + qJ(2);
t125 = cos(t127) / 0.2e1 + cos(t126) / 0.2e1;
t1 = [0, t136 -(-t133 * t125 - t134 * t132) * t128 + t130 * t136, 0, 0, 0; 0, -t135 -(t134 * t125 - t133 * t132) * t128 - t130 * t135, 0, 0, 0; 1, t131 -(sin(t126) / 0.2e1 + sin(t127) / 0.2e1) * t128 + t131 * t130, 0, 0, 0;];
Jg_rot  = t1;
