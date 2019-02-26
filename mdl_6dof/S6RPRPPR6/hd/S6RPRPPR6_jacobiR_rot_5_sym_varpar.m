% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPPR6_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:41:58
% EndTime: 2019-02-26 20:41:58
% DurationCPUTime: 0.05s
% Computational Cost: add. (25->11), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
t57 = sin(pkin(10));
t59 = sin(qJ(1));
t64 = t59 * t57;
t58 = cos(pkin(10));
t63 = t59 * t58;
t60 = cos(qJ(1));
t62 = t60 * t57;
t61 = t60 * t58;
t56 = qJ(3) + pkin(9);
t55 = cos(t56);
t54 = sin(t56);
t1 = [t54 * t61 - t64, 0, t55 * t63, 0, 0, 0; t54 * t63 + t62, 0, -t55 * t61, 0, 0, 0; 0, 0, -t54 * t58, 0, 0, 0; -t54 * t62 - t63, 0, -t55 * t64, 0, 0, 0; -t54 * t64 + t61, 0, t55 * t62, 0, 0, 0; 0, 0, t54 * t57, 0, 0, 0; -t60 * t55, 0, t59 * t54, 0, 0, 0; -t59 * t55, 0, -t60 * t54, 0, 0, 0; 0, 0, t55, 0, 0, 0;];
JR_rot  = t1;
