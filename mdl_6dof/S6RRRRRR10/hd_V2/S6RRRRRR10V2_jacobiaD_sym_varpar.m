% Zeitableitung der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6RRRRRR10V2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
% 
% Output:
% JaD [6x6]
%   Zeitableitung der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-06-09 18:09
% Revision: 5b4a470511447f972f018abf50a0dfc4f1a1159f (2019-06-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD = S6RRRRRR10V2_jacobiaD_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)



%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR10V2_jacobiaD_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR10V2_jacobiaD_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_sym_varpar: pkin has to be [6x1] (double)');
%% Function calls
if link_index == 0
	JaD=S6RRRRRR10V2_jacobiaD_0_sym_varpar(qJ, qJD, r_i_i_C, pkin);
elseif link_index == 1
	JaD=S6RRRRRR10V2_jacobiaD_1_sym_varpar(qJ, qJD, r_i_i_C, pkin);
elseif link_index == 2
	JaD=S6RRRRRR10V2_jacobiaD_2_sym_varpar(qJ, qJD, r_i_i_C, pkin);
elseif link_index == 3
	JaD=S6RRRRRR10V2_jacobiaD_3_sym_varpar(qJ, qJD, r_i_i_C, pkin);
elseif link_index == 4
	JaD=S6RRRRRR10V2_jacobiaD_4_sym_varpar(qJ, qJD, r_i_i_C, pkin);
elseif link_index == 5
	JaD=S6RRRRRR10V2_jacobiaD_5_sym_varpar(qJ, qJD, r_i_i_C, pkin);
elseif link_index == 6
	JaD=S6RRRRRR10V2_jacobiaD_6_sym_varpar(qJ, qJD, r_i_i_C, pkin);
else
	JaD=NaN(6,6);
end