% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
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
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR10_jacobiR_rot_2_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobiR_rot_2_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobiR_rot_2_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_2_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:17
% EndTime: 2018-11-23 11:27:17
% DurationCPUTime: 0.03s
% Computational Cost: add. (38->12), mult. (38->16), div. (0->0), fcn. (48->9), ass. (0->18)
t69 = pkin(6) - qJ(2);
t59 = pkin(6) + qJ(2);
t68 = sin(t59) / 0.2e1;
t67 = sin(t69);
t53 = t68 - t67 / 0.2e1;
t62 = sin(qJ(1));
t63 = cos(qJ(2));
t64 = cos(qJ(1));
t66 = -t64 * t53 - t62 * t63;
t65 = t62 * t53 - t64 * t63;
t61 = sin(qJ(2));
t60 = sin(pkin(6));
t58 = cos(t69);
t57 = cos(t59) / 0.2e1;
t54 = t58 / 0.2e1 + t57;
t52 = -t62 * t54 - t64 * t61;
t51 = -t64 * t54 + t62 * t61;
t1 = [t66, t52, 0, 0, 0, 0; -t65, -t51, 0, 0, 0, 0; 0, t68 + t67 / 0.2e1, 0, 0, 0, 0; t51, t65, 0, 0, 0, 0; t52, t66, 0, 0, 0, 0; 0, t57 - t58 / 0.2e1, 0, 0, 0, 0; t64 * t60, 0, 0, 0, 0, 0; t62 * t60, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JR_rot  = t1;
