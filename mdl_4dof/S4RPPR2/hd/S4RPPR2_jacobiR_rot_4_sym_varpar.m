% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4RPPR2
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
%
% Output:
% JR_rot [9x4]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S4RPPR2_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_jacobiR_rot_4_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_jacobiR_rot_4_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:31:22
% EndTime: 2019-02-26 19:31:22
% DurationCPUTime: 0.02s
% Computational Cost: add. (24->7), mult. (16->4), div. (0->0), fcn. (32->4), ass. (0->8)
t41 = cos(qJ(1));
t40 = sin(qJ(1));
t39 = pkin(6) + qJ(4);
t38 = cos(t39);
t37 = sin(t39);
t33 = t41 * t37 - t40 * t38;
t32 = -t40 * t37 - t41 * t38;
t1 = [t33, 0, 0, -t33; -t32, 0, 0, t32; 0, 0, 0, 0; -t32, 0, 0, t32; -t33, 0, 0, t33; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
JR_rot  = t1;
